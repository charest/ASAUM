#include "amr_flags.hpp"
#include "box_mesh.hpp"
#ifdef HAVE_EXODUS
#include "exodus_mesh.hpp"
#endif
#include "file_type.hpp"
#include "mesh_boundary.hpp"
#include "mesh_struct_node.hpp"
#include "mesh_struct_block.hpp"
#include "mesh_unstruct_block.hpp"
#include "partition.hpp"

#include "comm/comm_map.hpp"
#include "comm/comm_queue.hpp"
#include "comm/comm_utils.hpp"
#include "lua/lua_utils.hpp"
#include "utils/arguments.hpp"
#include "utils/cast.hpp"

namespace prl {

// command line options specific to the mesh
static option_category_t  OptionsCat("MESH OPTIONS:");

static option_t<std::string> MeshOption(
    "m",
    OptionsCat,
    option_description_t("Specify the mesh"),
    value_description_t("FILENAME"));

static option_t<std::string> DistOption(
    "d",
    OptionsCat,
    option_description_t(
      "Specify the file distribution method: "
      "sequential (default), cyclic, hostname."));

static option_t<std::string> PartOption(
    "a",
    OptionsCat,
    option_description_t(
      "Specify the partition algorithm: "
      "kway (default), geomkway, geom, metis, naive."));

  
////////////////////////////////////////////////////////////////////////////////
/// Check if there are anny changes
////////////////////////////////////////////////////////////////////////////////
bool has_changes(const std::vector<amr_flags_t> & refine_flags)
{
  bool has_change = false;
  for (auto & flg : refine_flags) {
    if (flg.has_change()) {
      has_change = true;
      break;
    }
  }
 return has_change;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor
////////////////////////////////////////////////////////////////////////////////
mesh_t::mesh_t(
    lua_result_t input,
    mpi_comm_t & comm,
    bool is_structured,
    const std::vector<int_t> & parts) :
  partitioning_(parts), comm_(comm) 
{

  auto is_root = comm.is_root();

  // read boundary conditions
  if (!input["boundaries"].empty()) {
    
    auto bnd_input = input["boundaries"];
    int_t num_bnds = bnd_input.size();
    
    boundaries_.reserve(num_bnds);
    
    for (int_t i=0; i<num_bnds; ++i) {
      
      auto this_bnd = bnd_input[i+1];
      auto res = make_mesh_boundary(this_bnd, comm_, i, is_structured);
      
      if (res) {
        const auto & name = res->name();  
        auto ret = boundary_map_.emplace(name, boundaries_.size());
        
        if (!ret.second) {
          if (is_root) {
            std::cout << "Duplicate boundary with name '" << name;
            std::cout << "' detected." << std::endl;
          }
          comm.exit(consts::failure);
        }

        boundaries_.emplace_back(std::move(res));
      }
    } // res

  }

}

////////////////////////////////////////////////////////////////////////////////
/// Get rank owner
////////////////////////////////////////////////////////////////////////////////
int_t mesh_t::block_owner(int_t bid) const
{ return owner( block_offsets_, bid); }

bool mesh_t::owns_block(int_t bid) const
{
  auto rank = comm_.rank();
  return bid >= block_offsets_[rank] && bid < block_offsets_[rank+1];
}

////////////////////////////////////////////////////////////////////////////////
/// Get the un/structured root
////////////////////////////////////////////////////////////////////////////////
mesh_struct_node_t * mesh_t::struct_root(int_t i)
{ return dynamic_cast<mesh_struct_node_t *>(roots_[i].get()); }

const mesh_struct_node_t * mesh_t::struct_root(int_t i) const
{ return dynamic_cast<const mesh_struct_node_t *>(roots_[i].get()); }

mesh_unstruct_node_t * mesh_t::unstruct_root(int_t i)
{ return dynamic_cast<mesh_unstruct_node_t *>(roots_[i].get()); }

const mesh_unstruct_node_t * mesh_t::unstruct_root(int_t i) const
{ return dynamic_cast<const mesh_unstruct_node_t *>(roots_[i].get()); }

////////////////////////////////////////////////////////////////////////////////
/// Get the un/structured block
////////////////////////////////////////////////////////////////////////////////
mesh_struct_block_t * mesh_t::struct_block(int_t i)
{ return dynamic_cast<mesh_struct_block_t *>(block(i)); }

const mesh_struct_block_t * mesh_t::struct_block(int_t i) const
{ return dynamic_cast<const mesh_struct_block_t *>(block(i)); }

mesh_unstruct_block_t * mesh_t::unstruct_block(int_t i)
{ return dynamic_cast<mesh_unstruct_block_t *>(block(i)); }

const mesh_unstruct_block_t * mesh_t::unstruct_block(int_t i) const
{ return dynamic_cast<const mesh_unstruct_block_t *>(block(i)); }


////////////////////////////////////////////////////////////////////////////////
/// Compute block counts
////////////////////////////////////////////////////////////////////////////////
std::vector<int_t> mesh_t::block_counts() const
{
  auto n = block_offsets_.size();
  if (!n) return {};

  std::vector<int_t> counts(n-1);
  for (size_t i=0; i<n-1; ++i)
    counts[i] = block_offsets_[i+1] - block_offsets_[i];

  return counts;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute block counts
////////////////////////////////////////////////////////////////////////////////
std::vector<int_t> mesh_t::block_ids() const
{
  std::vector<int_t> ids;
  ids.reserve(num_blocks());
  for (auto & b : blocks_)
    ids.emplace_back(b->id());
  return ids;
}

////////////////////////////////////////////////////////////////////////////////
/// Initialize a mesh
////////////////////////////////////////////////////////////////////////////////
void mesh_t::initialize()
{
  auto is_root = comm_.is_root();

  number();
  
  if (num_ghost_) {
    if (is_root) {
      std::cout << "msh> Expand mesh by " << num_ghost_ << " cells";
      std::cout << std::endl;
    }
    build_halo();
  }

  if (is_root) {
    std::cout << "msh> Initializing blocks";
    std::cout <<std::endl;
  }
  
  for (auto & b : blocks_)
    b->initialize(with_face_geom_, with_cell_geom_, boundaries_);
    
  if (num_ghost_ && with_cell_geom_) exchange_geometry();

}

////////////////////////////////////////////////////////////////////////////////
/// Build the halo
////////////////////////////////////////////////////////////////////////////////
void mesh_t::build_halo()
{
  comm_maps_.reset();
  if (num_ghost_) {
    std::vector<comm_map_block_t> maps;
    build_halo(maps);
    comm_maps_ = std::make_unique<comm_map_t>( std::move(maps) );
    ghost_exchange();
  }
}
  
////////////////////////////////////////////////////////////////////////////////
/// Exchange geometry
////////////////////////////////////////////////////////////////////////////////
void mesh_t::exchange_geometry()
{
  if (!(comm_maps_ && with_cell_geom_)) return;

  comm_queue_t q(comm_, *comm_maps_);
  int_t nblock = num_blocks();

  for (int_t b=0; b<nblock; ++b)
    blocks_[b]->queue_geometry_exchange(q[b]);

  for (int_t b=0; b<nblock; ++b)
    q[b].process();

  for (int_t b=0; b<nblock; ++b)
    q[b].wait();

}

 

////////////////////////////////////////////////////////////////////////////////
/// Make a mesh
////////////////////////////////////////////////////////////////////////////////
std::unique_ptr<mesh_t> make_mesh(
    lua_result_t inputs,
    mpi_comm_t & comm,
    const std::vector<int_t> & parts)
{
  if (inputs.empty()) {
    return {};
  }

  auto is_root = comm.is_root();
    
  // get mesh type
  auto valid = validate<std::string>(
      inputs,
      "type",
      {
        "structured",
#ifdef HAVE_EXODUS
       "exodus",
#endif
       "auto"
      },
      is_root);
  if (!valid) comm.exit(consts::failure);
  auto mesh_type = inputs["type"].as<std::string>();

  // mesh file specified on command line takes precedence
  std::string filename;
  if (MeshOption) {
    filename = MeshOption;
  }
  else if (!inputs["filename"].empty()) {
    auto res = as_scalar<std::string>(inputs, "filename", is_root);
    if (!res.exists) comm.exit(consts::failure);
    filename = res.value;
  }

  // determine the file distribution order
  distribution_t distribution = distribution_t::sequential;
  if (DistOption) {
    std::string str = DistOption;
    if (str == "sequential")
      distribution = distribution_t::sequential;
    else if (str == "cyclic")
      distribution = distribution_t::cyclic;
    else if (str == "hostname")
      distribution = distribution_t::hostname;
    else {
      if (is_root) {
        std::cout << "Unknown distribution option '" << str;
        std::cout << "'." << std::endl;
      }
      comm.exit(consts::failure);
    }
  }
  UNUSED(distribution);
  
  // determine the partition algorithm
  bool repart = false;
  part_alg_t part_alg = part_alg_t::parmetis_kway;
  if (PartOption) {
    std::string str = PartOption;
    if (str == "kway")
      part_alg = part_alg_t::parmetis_kway;
    else if (str == "geomkway")
      part_alg = part_alg_t::parmetis_geomkway;
    else if (str == "geom")
      part_alg = part_alg_t::parmetis_geom;
    else if (str == "metis")
      part_alg = part_alg_t::metis;
    else if (str == "naive")
      part_alg = part_alg_t::naive;
    else {
      if (is_root) {
        std::cout << "Unknown partition algorithm '" << str;
        std::cout << "'." << std::endl;
      }
      comm.exit(consts::failure);
    }
    repart = true;
  }
  UNUSED(repart);
  UNUSED(part_alg);
  
  
  // if mesh type is auto, try to figure it out.
  if  (mesh_type == "auto") {
    // need a filename to do the detection
    if (filename.empty()) {
      if (is_root) {
        std::cout << "Need a filename (either via input file or command line)";
        std::cout << " to do 'auto' detection." << std::endl;
      }
      comm.exit(consts::failure);
    }
    auto ft = get_file_type(filename);
    if (ft == file_type_t::partitioned_exodus ||
        ft == file_type_t::exodus) {
#ifdef HAVE_EXODUS
      mesh_type = "exodus";
#else
      if (is_root) {
        std::cout << "Exodus meshes not supported, rebuild with exodus";
        std::cout << " to support." << std::endl;
      }
      comm.exit(consts::failure);
#endif
    }
    else {
      if (is_root) {
        std::cout << "Could not ascertain the type of mesh that ";
        std::cout << " user provided file, '" << filename << "' is.";
        std::cout << "Try specifying 'type' to override.";
        std::cout << std::endl;
      }
      comm.exit(consts::failure);
    }
  }

  // load the determined one
  if (mesh_type == "structured") {
    return std::make_unique<box_mesh_t>(inputs, comm, parts);
  }
  else if (mesh_type == "exodus") {
#ifdef HAVE_EXODUS
    return std::make_unique<exodus_mesh_t>(
        inputs,
        comm,
        parts,
        filename,
        distribution,
        part_alg,
        repart);  
#else
    if (is_root) {
      std::cout << "Exodus meshes not supported, rebuild with exodus";
      std::cout << " to support." << std::endl;
    }
    comm.exit(consts::failure);
    return {};
#endif
  }

  return {};
}

////////////////////////////////////////////////////////////////////////////////
/// Find which cell a position is located in
////////////////////////////////////////////////////////////////////////////////
bool mesh_t::find_cell(const real_t * x, int_t & block_id, int_t & cell_id)
{

  auto nblock = num_blocks();

  for (block_id=0; block_id<nblock; ++block_id) {

    auto mesh_block = block(block_id);
    const auto & bbox = mesh_block->bounding_box_;

    // test bounding box first
    if (bounding_box_contains(bbox, x, num_dims_))
      // now search block
      if (mesh_block->find_cell(x, cell_id))
          return true;
        
  }

  return false;

}

////////////////////////////////////////////////////////////////////////////////
/// Find which cell a position is located in
////////////////////////////////////////////////////////////////////////////////
bool mesh_t::find_cell(long_t gid, int_t & block_id, int_t & cell_id)
{

  auto nblock = num_blocks();

  for (block_id=0; block_id<nblock; ++block_id) {

    if (blocks_[block_id]->find_cell(gid, cell_id))
      return true;
        
  }

  return false;

}


////////////////////////////////////////////////////////////////////////////////
/// Splice connectivity
////////////////////////////////////////////////////////////////////////////////
void mesh_t::ghost_exchange()
{

  int_t num_blks = comm_maps_->size();
  sub_comm_t subcomm(comm_, num_blks);
  if (!comm_maps_ || comm_maps_->empty()) return;
  
  auto world_rank = comm_.rank();
  auto comm_size = subcomm.comm().size();
  
  //----------------------------------------------------------------------------
  // Count
  //----------------------------------------------------------------------------

  // first count for storage
  std::vector<byte_t> tempbuf;
  
  std::vector<int> sendcounts(comm_size, 0);
  std::fill(sendcounts.begin(), sendcounts.end(), 0);

  std::map<int_t, int_t> block_map;

  //------------------------------------
  // For each block,rank combo, we will be sending:
  // 1. the source block id
  // 2. the destination block id
  // 3. the number of shared entities
  // 4. For each shared entity:
  //    ...  the number of connections, and their global ids
  for (int_t i=0; i<num_blks; ++i) {
    
    auto this_block = block(i);
    block_map[ this_block->id() ] = i;

    // get shared info
    const auto & cm = comm_maps_->operator[](i);
    const auto & shared_ids = cm.shared_ids_;
    const auto & shared_offsets = shared_ids.offsets;
    const auto & shared_indices = shared_ids.indices;
    const auto & shared_users = cm.shared_users_;
    int_t num_shared_blocks = shared_ids.size();

    // now count
    for (int_t b=0; b<num_shared_blocks; ++b) {
      auto r = shared_users[b];
      if (r != world_rank) {
        auto r2 = subcomm.map_rank(r);
        if (r2 < 0) THROW_ERROR("Rank maps to one with no data!");
        // shared id start/end
        auto start = shared_offsets[b];
        auto end = shared_offsets[b+1];
        // src/dst block ids
        sendcounts[r2] += 2*sizeof(int_t);
        // iterator over each shared id
        for (auto j=start; j<end; ++j) {
          auto id = shared_indices[j];
          tempbuf.clear();
          this_block->pack(id, tempbuf, true);
          // number of bytes plus buffer
          sendcounts[r2] += sizeof(int_t) + tempbuf.size();
        } // shared
      } // rank
    } // block 

  } // queues

  // finish displacements
  std::vector<int> senddispls(comm_size+1);
  senddispls[0] = 0;
  std::partial_sum(sendcounts.begin(), sendcounts.end(), &senddispls[1]);

  //----------------------------------------------------------------------------
  // now fill buffers
  //----------------------------------------------------------------------------
  
  std::vector<byte_t> sendbuf(senddispls[comm_size]);
  
  std::fill(sendcounts.begin(), sendcounts.end(), 0);
  
  for (int_t i=0; i<num_blks; ++i) {
    
    auto this_block = block(i);
    auto block_id = this_block->id();

    // get shared info
    const auto & cm = comm_maps_->operator[](i);
    const auto & shared_ids = cm.shared_ids_;
    const auto & shared_offsets = shared_ids.offsets;
    const auto & shared_indices = shared_ids.indices;
    const auto & shared_users = cm.shared_users_;
    const auto & shared_blocks = cm.shared_blocks_;
    int_t num_shared_blocks = shared_ids.size();

    // now count
    for (int_t b=0; b<num_shared_blocks; ++b) {
      auto r = shared_users[b];
      if (r != world_rank) {
        auto r2 = subcomm.map_rank(r);
        if (r2 < 0) THROW_ERROR("Rank maps to one with no data!");
    
        // get buffer offset
        auto buf = sendbuf.data() + senddispls[r2] + sendcounts[r2];
        auto buf_start = buf; 

        // add src/dst block
        cast(&block_id, 1, buf);
        cast(&shared_blocks[b], 1, buf);
        
        // Get begin/end of shared ids 
        auto start = shared_offsets[b];
        auto end = shared_offsets[b+1];
        
        // iterate over shared ids
        for (auto j=start; j<end; ++j) {
          // get data
          auto id = shared_indices[j];
          tempbuf.clear();
          this_block->pack(id, tempbuf, true);
          // copy into buffer
          int_t nbytes = tempbuf.size();
          cast(&nbytes, 1, buf);
          cast(tempbuf.data(), tempbuf.size(), buf);
        } // shared
        
        // increment counter
        sendcounts[r2] += buf - buf_start; 
     
      } // rank
    } // block 

  } // queues

  // verify counts
  for (decltype(comm_size) r=0; r<comm_size; ++r)
    if (senddispls[r]+sendcounts[r] != senddispls[r+1])
      THROW_ERROR("Send counts do not match the expected values.");

  //----------------------------------------------------------------------------
  // Send shared information
  //----------------------------------------------------------------------------
      
  // send the counts
  std::vector<int> recvcounts(comm_size);
  subcomm.comm().all_to_all(sendcounts, recvcounts);

  // how much info will we be receiving
  std::vector<int> recvdispls(comm_size + 1);
  recvdispls[0] = 0;
  std::partial_sum(recvcounts.begin(), recvcounts.end(), &recvdispls[1]);

  std::vector<byte_t> recvbuf(recvdispls[comm_size]);

  // now send the actual vertex info
  auto req = subcomm.comm().all_to_allv(
      sendbuf,
      sendcounts,
      senddispls,
      recvbuf,
      recvcounts,
      recvdispls);

  //----------------------------------------------------------------------------
  // Determine ordering
  //----------------------------------------------------------------------------

  //------------------------------------
  // Ghosts needing local data
  
  std::vector< std::vector< std::pair<int_t, int_t> > >
    ghost_positions(num_blks);

  for (int_t idst=0; idst<num_blks; ++idst) {

    auto dst_block = block(idst);
    auto num_cells = dst_block->num_owned_cells();
    
    const auto & dst_map  = comm_maps_->operator[](idst);
    auto dst_block_id = dst_map.block_id_;
    const auto & ghost_owners = dst_map.ghost_owners_;
    const auto & ghost_blocks = dst_map.ghost_blocks_;
    const auto & ghost_ids = dst_map.ghost_ids_;
    const auto & ghost_offsets = ghost_ids.offsets;
    const auto & ghost_indices = ghost_ids.indices;
    int_t num_blks = ghost_owners.size();
    int_t tot_ghost = ghost_indices.size();

    auto & ghost_pos = ghost_positions[idst];
    ghost_pos.resize(tot_ghost);
       
    for (int_t bdst=0; bdst<num_blks; ++bdst) {
      
      auto r = ghost_owners[bdst];
      if (r == world_rank) {
        
        auto src_block_id = ghost_blocks[bdst];
        auto isrc = block_map.at(src_block_id);
    
        const auto & src_map = comm_maps_->operator[](isrc);
        const auto & shared_ids = src_map.shared_ids_;
        const auto & shared_offsets = shared_ids.offsets;
        const auto & shared_indices = shared_ids.indices;

        auto bsrc = src_map.shared_block_map_.at(dst_block_id);

        auto shared_start = shared_offsets[bsrc];
        auto shared_end = shared_offsets[bsrc+1];
        auto num_shared = shared_end - shared_start;
      
        auto ghost_start = ghost_offsets[bdst];
        auto ghost_end = ghost_offsets[bdst+1];
        auto num_ghost = ghost_end - ghost_start;

        if (num_ghost != num_shared) {
          THROW_ERROR("Ghost and shared sizes do not match.  "
            << "block " << dst_block_id << " has " << num_ghost
            << " ghost but block " << src_block_id << " has "
            << num_shared);
        }

        for (int_t i=0; i<num_ghost; ++i) {
          auto src_id = shared_indices[shared_start + i];
          auto dst_id = ghost_indices[ghost_start + i];
          ghost_pos[ dst_id - num_cells ] = {isrc, src_id};
        } // i
      
      } // rank
    } // block
  
  } // maps
  
  //------------------------------------
  // Ghosts needing other ranks data
  
  // need data now
  req.wait();
  
  const auto * buf = recvbuf.data();
  const auto * buf_start = buf;
    
  for (int_t r=0; r<comm_size; ++r) { 
    auto buf_end = &recvbuf[ recvdispls[r+1] ];
    for (; buf<buf_end;) {

      // get source and destination block id
      int_t src_block_id;
      uncast(buf, 1, &src_block_id);
      
      int_t dst_block_id;
      uncast(buf, 1, &dst_block_id);
    
      // find the correct queue
      const auto & idst = block_map.at(dst_block_id);
      const auto & map  = comm_maps_->operator[](idst);
      const auto & ghost_owners = map.ghost_owners_;
      const auto & ghost_ids = map.ghost_ids_;
      const auto & ghost_offsets = ghost_ids.offsets;
      const auto & ghost_indices = ghost_ids.indices;
      int_t num_blks = ghost_owners.size();
    
      auto & ghost_pos = ghost_positions[idst];

      auto this_block = block(idst);
      auto num_cells = this_block->num_owned_cells();

      // and its destination block
      auto b = map.ghost_block_map_.at(src_block_id);
      if (b >= num_blks) {
        THROW_ERROR("The desitnation block is out of bounds.  "
            << "Block " << dst_block_id << " only has "
            << num_blks << " ghost blocks, you are asking for the "
            << b << "th block.");
      }
  
      // Unpack the data
      auto rnk = ghost_owners[b];
      if (rnk != r) {
        THROW_ERROR("The ghost block came from rank "
            << r << " but is supposed to have come from "
            << rnk << ".");
      }
        
      auto r2 = subcomm.map_rank(r);
      if (r2 < 0) THROW_ERROR("Rank maps to one with no data!");
        
      auto start = ghost_offsets[b];
      auto end = ghost_offsets[b+1];
    
      for (int_t i=start; i<end; ++i) {
        // get number of bytes
        int_t nbytes;
        uncast(buf, 1, &nbytes);
        // mark location
        auto id = ghost_indices[i];
        ghost_pos[ id - num_cells ] = {-1, buf - buf_start};
        // bump counter
        buf += nbytes;
      } // i

    } // displs

    if (buf != buf_end)
      THROW_ERROR("Receive counts do not match the expected values.");

  } // comm
  
  //----------------------------------------------------------------------------
  // Unpack data
  //----------------------------------------------------------------------------

  for (int_t idst=0; idst<num_blks; ++idst) {

    auto dst_block = block(idst);

    for (const auto & pair : ghost_positions[idst]) {
      auto isrc = pair.first;

      // Local exchange
      if (isrc > -1) {
        auto src_id = pair.second;
        auto src_block = block(isrc);
        tempbuf.clear();
        src_block->pack(src_id, tempbuf, true);
        const auto * buf = tempbuf.data();
        dst_block->unpack(buf, true);
      }
      // Global exchange
      else {
        auto offset = pair.second;
        const auto * buf = recvbuf.data() + offset;
        dst_block->unpack(buf, true);
      }
    
    } // ghosts

  } // blocks

}

////////////////////////////////////////////////////////////////////////////////
/// Compute global ids for disjoint inndexes
////////////////////////////////////////////////////////////////////////////////
void mesh_t::number_disjoint(int_t dim)
{

  auto nblocks = num_blocks();

  std::vector<long_t> offsets(nblocks+1);
  offsets[0] = 0;

  long_t cnt = 0;
  for (int_t b=0; b<nblocks; ++b) {
    auto n = blocks_[b]->num_owned_entities(dim);
    cnt += n;
    offsets[b+1] = offsets[b]+ n;
  }

  // construct the global mapping
  std::vector<long_t> dist;
  global_offsets(comm_, cnt, dist); 
  
  // now compute the indices
  auto comm_rank = comm_.rank();
  auto start = dist[comm_rank];

  for (int_t b=0; b<nblocks; ++b)
    blocks_[b]->number_disjoint(dim, start + offsets[b]);


}


////////////////////////////////////////////////////////////////////////////////
/// Number leaves
////////////////////////////////////////////////////////////////////////////////
void mesh_t::number_leaves() {
  int_t counter = 0;
  for (auto & b : roots_)
    b->number(counter);
}


////////////////////////////////////////////////////////////////////////////////
/// count leaves
////////////////////////////////////////////////////////////////////////////////
size_t mesh_t::count_leaves() const
{
  size_t n = 0;
  for (const auto & b : roots_)
    n += b->count_leaves();
  return n;
}

////////////////////////////////////////////////////////////////////////////////
/// extract leaves
////////////////////////////////////////////////////////////////////////////////
void mesh_t::get_leaves(std::vector<mesh_node_t*> & list) const
{
  auto n = count_leaves();
  list.clear();
  list.reserve(n);
  for (const auto & b : roots_)
    b->get_leaves(list);
}

std::vector<mesh_node_t*> mesh_t::get_leaves() const
{
  std::vector<mesh_node_t*> leaves;
  get_leaves(leaves);
  return leaves;
}

////////////////////////////////////////////////////////////////////////////////
/// extract blocks
////////////////////////////////////////////////////////////////////////////////
void mesh_t::get_leaf_blocks(std::vector<mesh_block_t*> & list) const
{
  auto n = count_leaves();
  list.clear();
  list.reserve(n);
  for (const auto & b : roots_)
    b->get_leaf_blocks(list);
}

std::vector<mesh_block_t*> mesh_t::get_leaf_blocks() const
{
  std::vector<mesh_block_t*> leaves;
  get_leaf_blocks(leaves);
  return leaves;
}
  

////////////////////////////////////////////////////////////////////////////////
// post amr setup
////////////////////////////////////////////////////////////////////////////////
void mesh_t::post_amr() {
  // re-mark as used/not used
  number_leaves();

  // reset block list
  get_leaf_blocks(blocks_);

  // update block offsets
  int_t nblocks = blocks_.size();
  global_offsets(comm_, nblocks, block_offsets_); 

  // now find new child neighbors
  find_neighbors();
}


////////////////////////////////////////////////////////////////////////////////
// Checck
////////////////////////////////////////////////////////////////////////////////
void mesh_t::fixup(
    std::vector<mesh_node_t*> & blocks,
    std::vector<amr_flags_t> & flags)
{
  bool first_pass = true;
  while (check(blocks, flags, first_pass))
    first_pass = false;

  change_connect(blocks, flags);
}

////////////////////////////////////////////////////////////////////////////////
// One check pass
////////////////////////////////////////////////////////////////////////////////
bool mesh_t::check(
    std::vector<mesh_node_t*> & blocks,
    std::vector<amr_flags_t> & flags,
    bool first_pass)
{
  
  bool is_changed = false;
  
  while (auto res =  inner_check(blocks, flags, first_pass)) {
    std::cout << "CHECK RES " << res << std::endl;
  }

  if (!first_pass) {
    auto ret = final_check(blocks, flags);
    if (ret) is_changed = true;
  }

  return first_pass ? true : is_changed;
}
  
////////////////////////////////////////////////////////////////////////////////
// One check pass
////////////////////////////////////////////////////////////////////////////////
bool mesh_t::inner_check(
  std::vector<mesh_node_t*> & blocks,
  std::vector<amr_flags_t> & flags,
  bool first_pass)
{
  bool is_changed = false;
  
  // Part 1 - Coarsening conflict check
  // 1. Sibling must have the same dimension in the direction of coarsening
  // 2. Sibling must also be flagged for coarsening.
  for (auto b : blocks) 
    is_changed = is_changed || b->first_sibling_check(&flags[0]);
  
  
  // Part 2 - Refinement level conflict check
  // 1. Eliminate inadmissible cases where the level difference > 2
  // 2. Set Nothing to refine or coarsen to nothnig in some cases
  for (auto b : blocks) 
    is_changed = is_changed || b->ratio_check(&flags[0], first_pass);
  
  // Part 3 - Max refinement level check
  for (auto b : blocks) 
    is_changed = is_changed || b->level_check(&flags[0], max_refinement_level_);

  return is_changed;
}

////////////////////////////////////////////////////////////////////////////////
// One check pass
////////////////////////////////////////////////////////////////////////////////
bool mesh_t::final_check(
    std::vector<mesh_node_t*> & blocks,
    std::vector<amr_flags_t> & flags)
{
  bool is_changed = false;

  // store old flags
  auto old_flags = flags;
  
  // Part 4 - Second sibling check
  // 1. Levels in the direction other than the coarsening direction
  //    are the same between block and sibling.

  // check siblings in direcctinos not alligned with direction of coarsening
  for (auto b : blocks) 
    is_changed = is_changed || b->final_check_1(&flags[0], &old_flags[0]);

  
  // Part 5 - Impossible cases check
  for (auto b : blocks) 
    is_changed = is_changed || b->final_check_2(&flags[0], &old_flags[0]);

  return is_changed;
}

////////////////////////////////////////////////////////////////////////////////
// Change connectivity
////////////////////////////////////////////////////////////////////////////////
void mesh_t::change_connect(
    std::vector<mesh_node_t*> & blocks,
    const std::vector<amr_flags_t> & flags)
{
  
  // rearrange connectivity if neihbor does not have a common parent
  for (auto b : blocks) 
    b->change_connect(&flags[0]);

}

////////////////////////////////////////////////////////////////////////////////
// Refine the mesh
////////////////////////////////////////////////////////////////////////////////
bool mesh_t::refine_tree(
    std::vector<mesh_node_t*> & nodes,
    const std::vector<amr_flags_t> & flags,
    amr_callback_t callback)
{
  bool did_refine = false;

  for (auto & n : nodes) {
    auto bid = n->id();
    const auto & rflags = flags[bid];

    auto res = n->refine_tree(rflags);
    did_refine = did_refine || res;

    auto block = n->block();

    if (res && block) {

      n->refine_block(rflags);
      
      if (callback) {
        auto new_blocks = n->get_leaf_blocks();
        callback(rflags, block, &new_blocks[0], new_blocks.size());
      }
      
      n->free_refined_blocks();

    } // res

  } // odes

  return did_refine;

}

////////////////////////////////////////////////////////////////////////////////
// Caorsen the mesh
////////////////////////////////////////////////////////////////////////////////
bool mesh_t::coarsen_tree(
    std::vector<mesh_node_t*> & nodes,
    const std::vector<amr_flags_t> & flags,
    amr_callback_t callback)
{
  bool did_coarsen = false;

  for (auto & n : nodes) {

    // cant coarsen roots
    if (!n->parent() || n->marked_for_deletion()) continue;

    auto bid = n->id();
    const auto & rflags = flags[bid];

    //----------------------------------
    // coarsen tree
    auto res = n->coarsen_tree(rflags);
    did_coarsen = did_coarsen || res;
    
    if (res) {

      auto block = n->block();
      auto parent = n->parent();
    
      if (block) {

        parent->coarsen_block();
        
        if (callback) {
          auto old_blocks = parent->get_leaf_blocks();
          assert(old_blocks.size() == 2);
          callback(rflags, parent->block(), &old_blocks[0], old_blocks.size());
        }
      
        parent->free_child_blocks();

      } // block

      parent->merge();
    } // res
    //----------------------------------

  } // blocks

  for (auto & rt : roots_)
    rt->free_marked_nodes();

  return did_coarsen;

}


////////////////////////////////////////////////////////////////////////////////
/// initial refinement of a mesh
////////////////////////////////////////////////////////////////////////////////
void mesh_t::initial_refinement(const  std::vector<int_t> & levels)
{

  //--- initialize the number of levels in each directoin
  std::vector<int_t> rlevels;
  int_t num_levels = levels.size();

  if (num_levels==1) {
    rlevels.assign(num_dims_, levels[0]);
  }
  else if (num_levels == num_dims_) {
    rlevels.assign(levels.begin(), levels.end());
  }
  else {
    if (comm_.is_root()) {
      std::cout << "Initial requires either one level for use in all directions, ";
      std::cout << "or " << num_dims_ << " levels (one for each dimension). ";
      std::cout << "You provided '" << num_levels << "'" << std::endl;
      std::cout << std::endl;
    }
    comm_.exit(consts::failure);
  }

  refine(rlevels.data());

  post_amr();

  number();
 
}

////////////////////////////////////////////////////////////////////////////////
// Main amr driver
////////////////////////////////////////////////////////////////////////////////
void mesh_t::refine(int_t * levels)
{
  bool refine_flags[] = {false, false, false};
  bool do_refine = false;

  for (int_t d=0; d<num_dims_; ++d) {
    refine_flags[d] = levels[d] > 0;
    do_refine = do_refine || refine_flags[d];
  }

  while (do_refine) {

    refine(refine_flags);

    do_refine = false;

    for (int_t d=0; d<num_dims_; ++d) {
      if (refine_flags[d])
        levels[d]--;
      refine_flags[d] = levels[d] > 0;
      do_refine = do_refine || refine_flags[d];
    }

  } // while

}

////////////////////////////////////////////////////////////////////////////////
// Main amr driver
////////////////////////////////////////////////////////////////////////////////
void mesh_t::refine(const bool * flags)
{
  
  // refine in all dims
  amr_flags_t refine_all;

  // refine in the specified dimensions
  if (flags) {
    for (int_t d=0; d<num_dims_; ++d)
      refine_all[d] = flags[d] ? amr_status_t::refine : amr_status_t::none;
  }
  // refine in all dimensions
  else  {
    for (int_t d=0; d<num_dims_; ++d)
      refine_all[d] = amr_status_t::refine;
  }


  if (refine_all.has_refine()) {

    // refine the tree
    auto nodes = get_leaves();
    for (auto & b : nodes) b->refine_tree(refine_all);
  
    // now refine my local blocks
    for (const auto n : nodes) {
      if (n->block()) {
        n->refine_block(refine_all);
        n->free_refined_blocks();
      }
    }

  } // has refine
  
}

////////////////////////////////////////////////////////////////////////////////
// Find new neighbors
////////////////////////////////////////////////////////////////////////////////
void mesh_t::find_neighbors()
{
  for (auto & rt : roots_)
    rt->clear_neighbors();
  
  for (auto & rt : roots_)
    rt->find_neighbors();
}

////////////////////////////////////////////////////////////////////////////////
// Main amr driver
////////////////////////////////////////////////////////////////////////////////
bool mesh_t::amr(
    const std::vector<amr_flags_t> & refine_flags,
    amr_callback_t refine_callback,
    amr_callback_t coarsen_callback)
{
  
  auto nodes = get_leaves();
  auto tot_blocks = nodes.size();

  //====================================
  // Gather flags

  auto comm_rank = comm_.rank();
  //auto comm_size = comm_.size();

  // gather block counts
  int_t num_blocks = refine_flags.size();
  auto block_cnts = block_counts();

  if (num_blocks != block_cnts[comm_rank]) 
    THROW_ERROR("Block counts and refine flags don't match!");

  // gather everyones refinement flags
  std::vector<amr_flags_t> global_refine_flags(tot_blocks);
  comm_.all_gatherv(refine_flags, global_refine_flags, block_cnts, block_offsets_);

  // make sure we have at least one change
  if (!has_changes(global_refine_flags)) return false;
 
  // Make sure refinement flags wont violate any rules
  std::cout << std::string(71, '=') << std::endl;
  for (size_t b=0; b<tot_blocks; ++b)
    std::cout << "block " << b << " marked " << global_refine_flags[b] 
      << " level " << nodes[b]->level(0) << std::endl;
  fixup(nodes, global_refine_flags);
  std::cout << std::endl;
  for (size_t b=0; b<tot_blocks; ++b)
    std::cout << "block " << b << " marked " << global_refine_flags[b] 
      << " level " << nodes[b]->level(0) << std::endl;
  std::cout << std::string(71, '=') << std::endl;
  
  // Check again after incase the fixup stopped any changes
  if (!has_changes(global_refine_flags)) return false;

  //====================================
  // Modify the tree
 
  auto did_refine = refine_tree(
      nodes,
      global_refine_flags,
      refine_callback);
  auto did_coarsen = coarsen_tree(
      nodes,
      global_refine_flags,
      coarsen_callback);
  
  if (!did_refine && !did_coarsen) return false;


  
  //====================================
  // finalize mesh blocks
  
  post_amr();

  number();
  
  build_halo();
  
  for (auto b :blocks_)
    b->initialize(with_face_geom_, with_cell_geom_, boundaries_);

  exchange_geometry();

  std::cout << "DID AMR!!!" << std::endl;

  return true;

}

} // namespace
