#include "mesh_builder.hpp"

#include "io/file_utils.hpp"
#include "io/outvar.hpp"
#include "mesh/amr_flags.hpp"
#include "mesh/mesh_struct_block.hpp"
#include "utils/arguments.hpp"
#include "utils/mark.hpp"
#include "utils/string_utils.hpp"
#include "utils/wall_time.hpp"

#include <numeric>
#include <set>

namespace prl {

// command line options specific to driver
option_category_t  OptionsCat("DRIVER OPTIONS:");

static option_t<std::string> CaseName(
    "c",
    OptionsCat,
    option_description_t("Case or output name"),
    value_description_t("CASENAME"));

static option_t<std::vector<int_t>> PartsOption(
    "p",
    OptionsCat,
    option_description_t("Specify number of partitions"),
    value_description_t("NPARTS"));

static option_t<std::vector<int_t>> RefineOption(
    "refine",
    OptionsCat,
    option_description_t("Specify number of initial refinement levels"),
    value_description_t("NLEVELS"));

////////////////////////////////////////////////////////////////////////////////
/// Constructor
////////////////////////////////////////////////////////////////////////////////
mesh_builder_t::mesh_builder_t(
    lua_t input,
    mpi_comm_t & comm) :
  input_(input), refine_(RefineOption), comm_(comm)
{
  auto is_root = comm_.is_root();

  // make sure there is a driver section
  if (!validate(input, "driver", is_root))
    comm.exit(consts::failure);
  auto driver_input = input["driver"];

  output_prefix_ = "out";
  if (CaseName) {
    output_prefix_ = CaseName;
  }
  else if (!driver_input["output_prefix"].empty()) {
    auto res = as_scalar<std::string>(driver_input, "output_prefix", is_root);
    if (!res.exists) comm.exit(consts::failure);
    output_prefix_ = res.value;
  }

  std::set<std::string> output_formats;
  
  if (!driver_input["output_formats"].empty()) {
    auto res = as_scalar<std::string>(driver_input, "output_formats", is_root);
    if (!res.exists) comm.exit(consts::failure);
    auto fmts = split(res.value, {' ', ',', ';'});
    output_formats.insert(fmts.begin(), fmts.end());
  }
    
  if (is_root && output_formats.size()) {
    std::cout << "drv> Output prefix: " << output_prefix_ << std::endl;
    std::cout << "drv> Output formats:";
    for (const auto & fmt : output_formats) std::cout << " " << fmt;
    std::cout << std::endl;
  }

  for (const auto & str : output_formats) {
    if (str != "vtk" && str != "csv" && str != "exo") {
      if (is_root)
        std::cout << "Output prefix of '" << str << "' is unsupported."  << std::endl;
      comm_.exit(consts::failure);
    }
  }

  auto vtk_input = input["vtk"];
  if (!vtk_input.empty()) {
    vtk_ = std::make_unique<vtk_writer_t>(vtk_input,  comm_);
  }
  else if (output_formats.count("vtk")) {
    vtk_ = std::make_unique<vtk_writer_t>();
  }
  
  auto csv_input = input["csv"];
  if (!csv_input.empty()) {
    csv_ = std::make_unique<csv_writer_t>(csv_input,  comm_);
  }
  else if (output_formats.count("csv")) {
    csv_ = std::make_unique<csv_writer_t>();
  }
  
  auto exo_input = input["exo"];
#ifdef HAVE_EXODUS
  if (!exo_input.empty()) {
    exo_ = std::make_unique<exo_writer_t>(exo_input,  comm_);
  }
  else if (output_formats.count("exo")) {
    exo_ = std::make_unique<exo_writer_t>();
  }
#else
  if (!exo_input.empty() || output_formats.count("exo")) {
    if (is_root) {
      std::cout << "Exodus output not supported, rebuild with exodus support ";
      std::cout << "to enable." << std::endl;
    }
    comm_.exit(consts::failure);
  }
#endif
  
    
  if (!validate(input_, "mesh", is_root)) comm_.exit(consts::failure);
  mesh_ = make_mesh(input["mesh"], comm_, PartsOption);
  if (!mesh_) comm_.exit(consts::failure);

  auto amr_input = input["amr"];
  if (!amr_input.empty()) {
    if (is_root) std::cout << "drv> AMR input found." << std::endl;
  }

  if (!driver_input["amr_levels"].empty()) {
    auto res = as_scalar<int_t>(driver_input, "amr_levels", is_root);
    if (!res.exists) comm.exit(consts::failure);
    amr_levels_ = res.value;
  }
  
  if (!driver_input["initial_amr"].empty()) {
    auto res = as_scalar<int_t>(driver_input, "initial_amr", is_root);
    if (!res.exists) comm.exit(consts::failure);
    initial_amr_ = res.value;
  }

  // setup comm stuff
  clock_start_ = wall_time();

}

////////////////////////////////////////////////////////////////////////////////
/// man driver unction
////////////////////////////////////////////////////////////////////////////////
void mesh_builder_t::build() {
  MARK_FUNCTION();
  
  auto is_root = comm_.is_root();

  // initial uniform refinement
  if (refine_.size()) {
    if (is_root) { 
      std::cout << "drv> Performing initial UNIFORM refinement: ";
      for (auto r : refine_) std::cout << r << " ";
      std::cout << std::endl;
    }
    mesh_->initial_refinement(refine_);
  }

  // initialize the mesh
  initialize();

  // amr
  mesh_->set_max_refinement_level(amr_levels_);
  if (initial_amr_>0) initial_amr();
  
  // output
  output();

}
  
////////////////////////////////////////////////////////////////////////////////
/// initialize the mesh
////////////////////////////////////////////////////////////////////////////////
void mesh_builder_t::initialize()
{
  MARK_FUNCTION();
  
  auto is_root = comm_.is_root();
  
  // create the requested data structures/connectivity/halo
  mesh_->build_halo();
  
  //-----------------------------------
  // build the mesh connectivity information
  auto num_dims = mesh_->num_dims();
  
  for (int_t i=0; i<mesh_->num_blocks(); ++i) {
    const auto mesh_block = mesh_->block(i);
    std::vector<std::pair<int_t, int_t>> conn;
      
    // cell to face
    conn.emplace_back(num_dims-1, num_dims);
    conn.emplace_back(num_dims, num_dims-1);
    
    // cell to vertex
    if (mesh_block->is_structured()) {
      conn.emplace_back(0, num_dims);
      conn.emplace_back(num_dims, 0);
    }
    // need face->vert for unstructured face geometry
    else if (num_dims>1) {
      conn.emplace_back(num_dims-1, 0);
      conn.emplace_back(num_dims, 0);
    }
  
    mesh_block->build_connectivity(conn);

    // build cell->cell neighbors through faces
    mesh_block->build_neighbors({std::make_pair(num_dims,num_dims-1)});

  } // blocks
  //-----------------------------------

  // prune any connectivity generated that wasnt requested
  mesh_->prune_connectivity();

  // compute geometric information
  mesh_->build_geometry(true, true);
  mesh_->exchange_geometry(); // for cell data

  // extract boundary information
  mesh_->build_boundaries();

}

  
////////////////////////////////////////////////////////////////////////////////
/// make a directory
////////////////////////////////////////////////////////////////////////////////
void mesh_builder_t::create_directory(const std::string & dirname)
{
  if (comm_.is_root()) {
    auto res = make_dir(dirname.c_str());
    
    if (res != 0 && res != EEXIST)
    {
      std::cout << "Unable to create directory '" << dirname << "', mkdir "
        << "returned errno=" << errno  << std::endl;
      comm_.exit(consts::failure);
    }
  }
  comm_.barrier();
}

////////////////////////////////////////////////////////////////////////////////
/// output the solution
////////////////////////////////////////////////////////////////////////////////
void mesh_builder_t::output()
{
  auto is_root = comm_.is_root();

  if (is_root) {
    std::cout << "drv> Output solution: " << output_prefix_ << " ("
      << output_counter_ << ")" << std::endl;
  }


  auto size = comm_.size();
  auto ndigits = num_digits(size);
  
  auto num_blocks = mesh_->num_blocks();
  auto tot_blocks = mesh_->tot_blocks();
  auto bdigits = num_digits(tot_blocks);

  auto iter_str = zero_padded(output_counter_);
  auto comm_size_str = zero_padded(tot_blocks, ndigits);
  auto block_size_str = zero_padded(tot_blocks, bdigits);

  std::stringstream title;
  title << "Solution snapshot at it=" << output_counter_
    << ", t=" << solution_time_;

  //----------------------------------------------------------------------------
  // VTK
  //----------------------------------------------------------------------------
  if (vtk_) {

    std::stringstream ss;
    
    ss << "vtkfiles";
    create_directory(ss.str());
    auto rootdir = ss.str();
    
    ss << "/" << iter_str;
    create_directory(ss.str());
    
    ss << "/" << output_prefix_;

    mesh_->output(
      vtk_.get(),
      ss.str(),
      output_counter_,
      solution_time_);

    std::stringstream vtm_ss, vtk_ss;
    vtm_ss << rootdir << "/" << output_prefix_;
    vtk_ss << iter_str << "/" << output_prefix_;
    
    mesh_->output_multiblock(
      vtk_.get(),
      vtm_ss.str(),
      vtk_ss.str(),
      output_counter_,
      solution_time_);

  } // vtk
  
  //----------------------------------------------------------------------------
  // CSV
  //----------------------------------------------------------------------------
  if (csv_) {

    std::stringstream ss;

    ss << "csvfiles";
    create_directory(ss.str());
    
    ss << "/" << iter_str;
    create_directory(ss.str());

    ss << "/" << output_prefix_ << ".mesh";
    mesh_->output(csv_.get(), ss.str(), title.str());

  } // csv
  
  //----------------------------------------------------------------------------
  // Exodus
  //----------------------------------------------------------------------------
#ifdef HAVE_EXODUS
  if (exo_) {

    std::stringstream ss;
    ss << "exofiles";
    create_directory(ss.str());

    ss << "/" << output_prefix_;
    
    mesh_->output(
      exo_.get(),
      ss.str(),
      title.str(),
      output_counter_,
      solution_time_);

  } // exo
#endif
  //----------------------------------------------------------------------------


  output_counter_++;

}

////////////////////////////////////////////////////////////////////////////////
/// Perform amr
////////////////////////////////////////////////////////////////////////////////
void mesh_builder_t::amr()
{
  
  auto is_root = comm_.is_root();
  if (is_root) std::cout << "drv> Performing amr." << std::endl;
  
  auto num_blocks = mesh_->num_blocks();
  
  std::vector<amr_flags_t> flags(num_blocks);
  //hydro_->flag_blocks_for_amr(flags.data());

  auto did_amr = mesh_->amr(
      flags,
      //--- Refine callback
      [&,this](auto flags, auto parent, auto children, auto n)
      { refine(flags, parent, children, n); },
      //--- Coarsen callback
      [&,this](auto flags, auto parent, auto children, auto n)
      { coarsen(flags, parent, children, n); }
  );

}

////////////////////////////////////////////////////////////////////////////////
/// Perform amr
////////////////////////////////////////////////////////////////////////////////
void mesh_builder_t::initial_amr()
{
  for (int_t i=0; i<initial_amr_; ++i)
    amr();
}

////////////////////////////////////////////////////////////////////////////////
/// Perform refinement
////////////////////////////////////////////////////////////////////////////////
void mesh_builder_t::refine(
    const amr_flags_t & flags,
    mesh_block_t * parent,
    mesh_block_t ** children,
    int_t nchild)
{
          
  //hydro_->refine(flags, parent, children, nchild);
  
  // clear memory
  auto bid = mesh_->block_global2local(parent->id());
}
  
////////////////////////////////////////////////////////////////////////////////
/// Perform Coarsening
////////////////////////////////////////////////////////////////////////////////
void mesh_builder_t::coarsen(
    const amr_flags_t & flags,
    mesh_block_t * parent,
    mesh_block_t ** children,
    int_t nchild)
{
  //hydro_->coarsen(flags, parent, children, nchild);
  std::cout << "NCCHILLD " << nchild << std::endl;
  
}
  
  
} // namespace

