#include "box_mesh.hpp"
#include "mesh_struct_node.hpp"
#include "mesh_struct_block.hpp"

#include "comm/comm_map.hpp"
#include "comm/comm_queue.hpp"
#include "comm/comm_utils.hpp"
#include "comm/entity_info.hpp"
#include "math/morton.hpp"
#include "math/decomp.hpp"
#include "math/subdivide.hpp"
#include "lua/lua_utils.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Box mesh constructor
////////////////////////////////////////////////////////////////////////////////
box_mesh_t::box_mesh_t(
    lua_result_t input,
    mpi_comm_t & comm,
    const std::vector<int_t> & parts) :
  mesh_t(input, comm, true, parts)
{
  // get communicator info
  auto is_root = comm.is_root();

  //----------------------------------------------------------------------------
  // Block coordinates
  
  // get block coordinates
  if (!validate(input, "coordinates", is_root))
    comm.exit(consts::failure);

  auto coords_input = input["coordinates"];
  auto ncoords = coords_input.size();
  if (ncoords<2) {
    if (is_root)
      std::cout << "Need at least two coordinates to form a block." << std::endl;
    comm.exit(consts::failure);
  }

  num_dims_ = 0;
  std::vector<real_t> block_vertices; 

  for (decltype(ncoords) i=0; i<ncoords; ++i) {
    auto inp = coords_input[i+1];
    auto coords = inp.as<std::vector<real_t>>();
    int_t n = coords.size();
    // invalid coordinates
    if (n < 1) {
      if (is_root)
        std::cout << "Invalid or zero-length block coordinates provided." << std::endl;
      comm.exit(consts::failure);
    }
    // if this is first time through, resize as necessary
    if (!num_dims_) {
      num_dims_ = n;
    }
    // make sure dims match
    else if (n != num_dims_) {
      if (is_root) {
        std::cout << "Each set of block coordinates must have the same dimensionality.";
        std::cout << "  Detected dimension was " << num_dims_;
        std::cout << ", but found coordinates of length " << n << std::endl;
      }
      comm.exit(consts::failure);
    }

    // add to list
    block_vertices.insert(block_vertices.end(), coords.begin(), coords.end());
  } // coords
 
    
  if (num_dims_ < 1 || num_dims_ > 3) {
    if (is_root) {
      std::cout << "msh> Mesh must be 1-, 2-, or 3-dimensional, you asked for a "
        << num_dims_ << "-dimensional mesh." << std::endl;
    }
    comm.exit(consts::failure);
  }
    
  int_t nparts = parts.size();
  if (!parts.empty() && nparts != 1) {
    if (is_root) {
      std::cout << "Structured meshes can only be partitioned by one dimension.";
      std::cout << "You provided " << nparts << " values.  Only one is allowed.";
      std::cout << std::endl;
    }
    comm.exit(consts::failure);
  }
  
  //----------------------------------------------------------------------------
  // Transformations

  //--- One provided
  if (!input["rotate"].empty()) {
  
    // get transformation
    auto ctm = as_vector<int_t>(input, "rotate", is_root);
    if (ctm.empty()) comm.exit(consts::failure);
    if (ctm.size() != static_cast<size_t>(num_dims_)) {
      if (is_root) {
        std::cout << "Provided otation matrix must be of length " << num_dims_;
        std::cout << ", you provided one of length " << ctm.size() << std::endl;
      }
      comm.exit(consts::failure);
    }

    // check bounds
    std::set<size_t> entries;
    for (auto i : ctm) {
      auto ii = std::abs(i);
      if (ii < 1 || ii > num_dims_+1) {
        if (is_root) {
          std::cout << "Absolute value of rotation entries must be between 1 and ";
          std::cout << num_dims_+1 << std::endl;
        }
        comm.exit(consts::failure);
      }
      if (!entries.emplace(ii).second) {
        if (is_root) {
          std::cout << "Duplicate entries in rotation detected." << std::endl;
        }
        comm.exit(consts::failure);
      }
    }

    // set transform
    transform_t tm(ctm.data(), num_dims_);

  
    // transform coords
    for (decltype(ncoords) v=0; v<ncoords; ++v)
      tm.transform( &block_vertices[v*num_dims_] );

  }
  
  //----------------------------------------------------------------------------
  // Initial mesh
  
  // create the initial mesh
  initial_mesh_ = std::make_unique<mesh_ijk_t>(num_dims_);
  initial_mesh_->reserve_vertices(ncoords);

  // add the coordinates
  for (decltype(ncoords) v=0; v<ncoords; ++v)
    initial_mesh_->add_vertex(&block_vertices[v*num_dims_]);
  
  //----------------------------------------------------------------------------
  // Blocks
  
  // get list of blocks
  if (!validate(input, "blocks", is_root))
    comm.exit(consts::failure);

  auto blocks_input = input["blocks"];
  auto nblocks = blocks_input.size();
  if (!nblocks) {
    if (is_root)
      std::cout << "Need at least one block." << std::endl;
    comm.exit(consts::failure);
  }

  initial_mesh_->reserve_blocks(nblocks);

  for (decltype(nblocks) b=0; b<nblocks; ++b) {
    auto block_input = blocks_input[b+1];

    //-----------------------------------
    // get dimensional info
    auto dims = as_vector<long_t>(block_input, "dimensions", is_root);
    if (dims.empty()) comm.exit(consts::failure);
    
    if (is_root) {
      std::cout << "msh> Creating block of size: ";
      for (auto d : dims) std::cout << d << " ";
      std::cout << "cells." << std::endl;
    }
    
    //-----------------------------------
    // get connectivity
    auto verts = as_vector<int_t>(block_input, "connect", is_root);
    if (verts.empty()) comm.exit(consts::failure);

    //--- one vertex per face
    if (num_dims_ == 1) {
      if (verts.size() != 2) {
        if (is_root) {
          std::cout << "Block connectivity in one-dimensions must be of size 2.";
          std::cout << "  You provided a list of length " << verts.size();
          std::cout << std::endl;
        }
        comm.exit(consts::failure);
      }
    }
    //--- two vertex per face
    else if (num_dims_ == 2) {
      if (verts.size() != 4) {
        if (is_root) {
          std::cout << "Block connectivity in two-dimensions must be of size 4.";
          std::cout << "  You provided a list of length " << verts.size();
          std::cout << std::endl;
        }
        comm.exit(consts::failure);
      }
    }
    //--- four vertex per face
    else {
      if (verts.size() != 8) {
        if (is_root) {
          std::cout << "Block connectivity in three-dimensions must be of size 8.";
          std::cout << "  You provided a list of length " << verts.size();
          std::cout << std::endl;
        }
        comm.exit(consts::failure);
      }
    }
   
    //-----------------------------------
    // Create block 
    
    auto & this_block = initial_mesh_->add_block(&dims[0], &verts[0], verts.size());
  
    //-----------------------------------
    // determine spacing function

    std::string spacing_type = "linear";

    if (!block_input["spacing"].empty()) {
      if (!validate<std::string>(block_input, "spacing", {"linear"}, is_root))
        comm.exit(consts::failure);

      spacing_type = block_input["spacing"].as<std::string>();
    }
      
    if (spacing_type == "linear")
      this_block.set_linear_spacing();
 
 
  } // blocks

  initial_mesh_->compute_neighbors();
  
}


////////////////////////////////////////////////////////////////////////////////
/// Load the mesh
////////////////////////////////////////////////////////////////////////////////
void box_mesh_t::load()
{
  auto is_root = comm_.is_root();
  auto comm_size = comm_.size();
  auto comm_rank = comm_.rank();
  
  if (is_root) std::cout << "msh> Generating mesh..." << std::endl;
  
  //----------------------------------------------------------------------------
  // partition mesh
  
  // determine how many blocks we have
  int_t tot_blocks = initial_mesh_->count_leaves();

  // determine desired partitioniing
  int_t nparts = partitioning_.empty() ? comm_size : partitioning_[0];
  
  // if parttionnig  doesnt  match (note, cant shrink number partitions)
  if (nparts > tot_blocks) {

    if (is_root) {
      std::cout << "msh> Splitting into ";
      std::cout << nparts << " blocks.";
      std::cout << std::endl;
    }

    if (initial_mesh_->partition(nparts)) {
      if (is_root) {
        std::cout << "Trouble partitioning! Could not acheive desired number";
        std::cout << " of partitions." << std::endl;
      }
      comm_.exit(consts::failure);
    }
  
  }
 
  // get the partitioned mesh
  auto initial_blocks = initial_mesh_->get_leaves();
  tot_blocks = initial_blocks.size();
  
  //----------------------------------------------------------------------------
  // Extract global mesh
  
  // create the roots
  reserve_roots(tot_blocks);
  for (int_t i=0; i<tot_blocks; ++i)
    add_root( std::make_unique<mesh_struct_node_t>(
          initial_blocks[i]->dims(), num_dims_, i) );

  // add neighbor info
  for (int_t i=0; i<tot_blocks; ++i) {
    
    for (const auto & sector_neigh : initial_blocks[i]->neighbors()) {
      const auto & sector = sector_neigh.first;
      for (const auto & neigh : sector_neigh.second)
        if (neigh.block().is_active()) {
          struct_root(i)->add_neighbor(
              sector,
              root(neigh.block().id()),
              neigh.ctm_me2donor_cells(),
              neigh.ctm_me2donor_verts(),
              neigh.range());
          struct_root(i)->add_root_neighbor(
              sector,
              root(neigh.block().id()),
              neigh.ctm_me2donor_cells(),
              neigh.ctm_me2donor_verts(),
              neigh.range());
        }
    }
    
    for (const auto & neigh : initial_blocks[i]->shared()) {
      if (neigh.block().is_active())
        struct_root(i)->add_shared(
            root(neigh.block().id()),
            neigh.sector(),
            neigh.range());
    }

  }
  
  //----------------------------------------------------------------------------
  // Distribute blocks

  // determine how many blocks we have
  if (tot_blocks < comm_size) {
    if (is_root) std::cout << "msh> Note: some ranks will not have blocks." << std::endl;
  }
  else if (tot_blocks > comm_size) {
    if (is_root) std::cout << "msh> Note: some ranks will have multiple blocks." << std::endl;
  }
  
  // subdivide blocks
  subdivide(tot_blocks, comm_size, block_offsets_);
  auto block_start = block_offsets_[comm_rank];
  auto num_blocks = block_offsets_[comm_rank+1] - block_start;
  
  //----------------------------------------------------------------------------
  // build local blocks
  
  for (int_t b=0; b<num_blocks; ++b) {

    // create the block
    auto bid = block_start + b;
    const auto & block_ijk = *initial_blocks[bid];
 
    auto this_block = std::make_unique<mesh_struct_block_t>(block_ijk);
    
    auto err = this_block->install_boundaries(block_ijk, boundaries_);
    if (err) comm_.exit(consts::failure);

    
    blocks_.emplace_back(this_block.get());
    root(bid)->set_block( std::move(this_block) );

  } // blocks
 
  // done with initial mesh
  initial_mesh_.reset();

}

////////////////////////////////////////////////////////////////////////////////
/// Construct the global mapping
////////////////////////////////////////////////////////////////////////////////
void box_mesh_t::number()
{
  number_disjoint(num_dims_);
  number_vertices();
}

////////////////////////////////////////////////////////////////////////////////
/// Add a halo
////////////////////////////////////////////////////////////////////////////////
void box_mesh_t::build_halo(std::vector<comm_map_block_t> & comm_maps)
{
  if (num_ghost_<0) return;
  
  //----------------------------------------------------------------------------
  // Determine cell distribution

  // get local and global number of blocks
  auto num_blks = num_blocks();
  auto tot_blks = tot_blocks();

  // determine my block sizes
  std::vector<long_t> block_sizes;
  block_sizes.reserve(num_blks);
  for (const auto & b : blocks_)
    block_sizes.emplace_back(b->num_owned_cells());

  // determine recv counts
  auto recvcounts = block_counts();

  // compute final distribution
  std::vector<long_t> cell_dist(tot_blks+1);
  global_offsets(comm_, block_sizes, cell_dist, recvcounts, block_offsets_);

  //----------------------------------------------------------------------------
  // Determine ghost / shared

  // resize the lists
  comm_maps.clear();
  comm_maps.reserve(num_blks);
    
  std::vector<entity_info_t> ghost_info;
  std::vector<entity_info_t> shared_info;
  
  for (int_t b=0; b<num_blks; ++b) {
    
    shared_info.clear();
    ghost_info.clear();

    auto this_block = dynamic_cast<mesh_struct_block_t*>(blocks_[b]);
    this_block->build_halo(
        num_ghost_,
        with_corners_,
        block_offsets_,
        cell_dist,
        shared_info,
        ghost_info); 

    // create a new comm map
    comm_maps.emplace_back(
        this_block->id(),
        this_block->num_all_cells(),
        shared_info,
        ghost_info);
    
  } // blocks

}
   
////////////////////////////////////////////////////////////////////////////////
/// Number the vertices
////////////////////////////////////////////////////////////////////////////////
void box_mesh_t::number_vertices()
{
  auto nblocks = num_blocks();
  
  std::vector<long_t> vert_counts(nblocks);

  for (int_t b=0; b<nblocks; ++b) {
    auto this_block = dynamic_cast<mesh_struct_block_t*>(blocks_[b]);
    auto unique_verts = this_block->count_vertices();
    vert_counts[b] = unique_verts;
  }
    
  
  // determine recv counts
  auto recvcounts = block_counts();

  // Vertices
  auto nblocks_tot = tot_blocks();
  std::vector<long_t> vert_dist(nblocks_tot+1);
  global_offsets(comm_, vert_counts, vert_dist, recvcounts, block_offsets_);
  
  //------------------------------------
  // Loop over blocks 
 
  std::vector<comm_map_block_t> comm_maps;
  comm_maps.reserve(nblocks);
    
  std::vector<entity_info_t> ghost_info;
  std::vector<entity_info_t> shared_info;
  
  for (auto & b : blocks_) {
    auto this_block = dynamic_cast<mesh_struct_block_t*>(b);

    // number the vertices that we own
    shared_info.clear();
    ghost_info.clear();
    this_block->number_vertices(
        block_offsets_,
        vert_dist[b->id()],
        shared_info,
        ghost_info);

    // create a new comm map
    comm_maps.emplace_back(
        this_block->id(),
        this_block->global_ids(0).size(),
        shared_info,
        ghost_info);

  } // blocks

  // create the map and exchange the vertex ids
  comm_map_t comm_map(std::move(comm_maps));
  comm_queue_t queue(comm_, comm_map);
  
  auto qit = queue.begin();
  for (auto & b : blocks_) {
    auto this_block = dynamic_cast<mesh_struct_block_t*>(b);
    qit->add(this_block->global_ids(0));
    qit++;
  }

  queue.process();
}

} // namespace
