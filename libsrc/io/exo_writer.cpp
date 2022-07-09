// user includes
#include "exo_writer.hpp"
#include "outvar.hpp"
#include "comm/mpi_comm.hpp"
#include "lua/lua_utils.hpp"

namespace prl {

//==============================================================================
//! \brief Constructor with input
//==============================================================================
exo_writer_t::exo_writer_t(lua_result_t input, mpi_comm_t & comm) 
{
  auto is_root = comm.is_root();
  
  // single or double figures
  if (!input["precision"].empty()) {
    
    auto res = validate<std::string>(
        input,
        "precision",
        {"single", "double"},
      is_root);
    if (!res) comm.exit(consts::failure);
    auto str = input["precision"].as<std::string>();

    if (str == "single")
      precision_ = precision_t::sp;
    else
      precision_ = precision_t::dp;
  }
}

//==============================================================================
//! \brief Open a file
//==============================================================================
int exo_writer_t::open(const std::string & name)
{
  
  // useful for debug
  ex_opts(EX_ABORT | EX_VERBOSE);
    
  // size of floating point variables used in app.
  int app_word_size = sizeof(real_t);
      
  // size of floating point to be stored in file.
  // change to float to save space
  int exo_word_size = sizeof(real_t);
  if (precision_ == precision_t::sp)
    app_word_size = sizeof(float);
  else if (precision_ == precision_t::dp)
    app_word_size = sizeof(double);

  // determine the file creation mode
  int cmode = EX_CLOBBER;

  // create file
  exoid_ = ex_create(name.c_str(), cmode, &app_word_size, &exo_word_size);
  if (exoid_ < 0) return consts::failure;

  return consts::success;
}
  
//==============================================================================
//! \brief Close a file
//==============================================================================
int exo_writer_t::close( void ) 
{ return ex_close(exoid_); }

////////////////////////////////////////////////////////////////////////////////
/// check if this is a large integer support
////////////////////////////////////////////////////////////////////////////////
bool exo_writer_t::is_int64()
{ return (ex_int64_status(exoid_) & EX_IDS_INT64_API); }
  
//============================================================================
//! \brief Helper function to make and initialize a set of exodus parameters.
//! \return the exodus parameters
//============================================================================
ex_init_params exo_writer_t::make_params(const char * label) {
  ex_init_params exopar;
  std::memset(&exopar, 0, sizeof(ex_init_params));
  if (label) std::strcpy(exopar.title, label);
  return exopar;
}

//============================================================================
//! \brief write the exodus parameters from a file.
//! \param [in] exo_id  The exodus file id.
//! \param [in] exo_params  The  exodus parameters
//============================================================================
int exo_writer_t::write_params(const ex_init_params & exo_params) {


  // keep track of the number of dimensions
  num_dims_ = exo_params.num_dim;

  // put the initialization parameters
  auto status = ex_put_init_ext(exoid_, &exo_params);
  if (status) return status;
  
  return 0;
}
    
//============================================================================
//! \brief Output variable names.
//============================================================================
int exo_writer_t::write_field_names(const std::vector<outvar_t> & vars)
{

  // count the variables/components
  int_t num_vars = vars.size();
  int_t var_counter = 0;

  var_map_.clear();
  for(int_t i=0; i<num_vars; ++i) {
    auto ncomp = vars[i].components;
    var_map_[vars[i].name] = varinfo_t{var_counter, ncomp};
    var_counter += ncomp;
  }

  // put the number of element fields
  auto status = ex_put_var_param(exoid_, "e", var_counter);
  if (status) return status;

  const char * var_ext[] = { "x", "y", "z" };

  // now register their names
  var_counter =  0;

  for(int_t i=0; i<num_vars; ++i)
  {
    auto ncomp = vars[i].components;
    auto dimensional = vars[i].dimensional;

    // multiple comonents
    if (ncomp>1) {
      for (int_t d=0; d<ncomp; ++d) {
        std::stringstream ss;
        ss << vars[i].name;
        if (dimensional) ss << var_ext[d];
        else ss << d;
        status = ex_put_var_name(exoid_, "e", var_counter + d+1, ss.str().c_str());
        if (status) return status;
      }
    }
    // scalar vaue
    else {
      status = ex_put_var_name(exoid_, "e", var_counter+1, vars[i].name.c_str());
      if (status) return status;
    }

    var_counter += ncomp;
  } // vars

  return 0;
}
        
//============================================================================
//! \brief Output variable
//============================================================================
int exo_writer_t::write_field(
    const real_t * data,
    size_t len,
    const std::string & name,
    int_t blk_id)
{

  auto it = var_map_.find(name);
  if (it == var_map_.end()) return -1;

  int time_step = 1;

  auto ncomp = it->second.components;
  auto id = it->second.id+1;

  // vectors
  if (ncomp > 1) {
    
    auto nvals = len / ncomp;
    std::vector<real_t> vals(nvals);

    for (int_t d=0; d<ncomp; ++d) {
      // copy data
      for (size_t i=0; i<nvals; ++i)
        vals[i] = data[i*ncomp + d];
      // plot
      auto status = ex_put_elem_var(
          exoid_,
          time_step,
          id + d,
          blk_id,
          nvals,
          vals.data());
      if (status) return status;
    }

  }
  // scalars
  else {
    auto status = ex_put_elem_var(
        exoid_,
        time_step,
        id,
        blk_id,
        len,
        data);
    if (status) return status;
  } // components
  
  return 0;
}

//============================================================================
//! \brief write the time
//============================================================================
int exo_writer_t::write_time(size_t time_step, real_t time)
{
  auto status = ex_put_time(exoid_, time_step, &time);
  if (status) return status;
  return 0;
}
  
//============================================================================
//! Write side set param only
//============================================================================
int exo_writer_t::write_side_set_param(
    int_t ss_id,
    size_t num_sides)
{
  auto status = ex_put_side_set_param(exoid_, ss_id, num_sides, 0);
  if (status) return status;
  return 0;
}

//============================================================================
//! \brief Output side sett names.
//============================================================================
int exo_writer_t::write_side_set_names(const std::vector<std::string> & names)
{
  auto num_names = names.size();
  std::vector<char *> names_cstr(num_names);
  for (size_t i=0; i<num_names; ++i) 
    names_cstr[i] = const_cast<char*>(names[i].c_str());
  
  auto status = ex_put_names(exoid_, EX_SIDE_SET, names_cstr.data());
  if (status) return status;
  
  return 0;
}
        

} // namespace prl
