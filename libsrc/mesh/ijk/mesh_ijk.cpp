#include "mesh_ijk.hpp"

#include "transform.hpp"
#include "utils/errors.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>

namespace prl {


////////////////////////////////////////////////////////////////////////////////
/// create a new mesh from an existing one 
////////////////////////////////////////////////////////////////////////////////
bool mesh_ijk_t::partition(int_t nparts, int_t min_size)
{

  const auto & permutations = sector_permutations[num_dims_-1];

  // get original leaves
  auto orig_blocks = get_leaves();
  int_t nblock = orig_blocks.size();

  //----------------------------------------------------------------------------
  // partition

  // parition blocks
  while (count_leaves() < static_cast<size_t>(nparts))
  {
    auto b = const_cast<block_ijk_t*>(find_largest_block());
    if (!b) return true;
    if (b->split_block(min_size)) return true;
  }
  
  // number the tree blocks
  number();

  //----------------------------------------------------------------------------
  // Find internal neighbors
  
  for (int_t b=0; b<nblock; ++b) {

    const auto orig_block = orig_blocks[b];
    
    // get new leaves
    auto new_blocks = orig_block->get_leaves();
    auto num_new_blocks = new_blocks.size();

    for (size_t b=0; b<num_new_blocks; ++b) {

      auto new_block = const_cast<block_ijk_t*>(new_blocks[b]);

      for (const auto & sector : permutations) {

        auto donor_blocks = new_block->find_neighbors(sector);
        
        auto target_range_copy = new_block->range_wrt_root(sector);
        
        for (auto d : donor_blocks) {
          auto donor_block = const_cast<block_ijk_t*>(d);
          auto target_range = target_range_copy;
          auto donor_range = donor_block->range();
          auto res = intersect_ranges(
              num_dims_,
              sector,
              target_range,
              donor_range);
          if (res) {
            new_block->add_internal_neighbor(
                *donor_block,
                sector,
                target_range);
          }
        } // donors
      } // permutations
    } // new blocks

  } // orig blocks
  
  //----------------------------------------------------------------------------
  // Find external neighbors
  for (int_t b=0; b<nblock; ++b) {

    const auto orig_block = orig_blocks[b];
    
    // copy cause neighbors might be changing
    auto neighbors = orig_block->neighbors();

    // Add out-of-root neighbors
    for (const auto & neigh_pair : neighbors) {
      
      const auto & target_sector = neigh_pair.first;
      const auto & sector_neighbors = neigh_pair.second;
        
      auto target_list = orig_block->find_sectors(target_sector);

      for (const auto & orig_neigh : sector_neighbors) {
        const auto & donor_node = orig_neigh.block();

        // want donor->me transformation
        auto target2donor = orig_neigh.tm_me2donor_cells();
        auto donor2target = reverse(target2donor);
              
        // transform ranges, looking for intersectoins 
        for (auto t : target_list) {
          auto target_node = const_cast<block_ijk_t*>(t);

          for (const auto & target_search_sector : permutations) {
            
            // skip irrelevant sectors
            if (!target_sector.derives_from(target_search_sector)) continue;
            
            // new target sector
            auto target_range = target_node->range_wrt_root(target_search_sector);

            // trasform it to donor frame
            auto target_range_in_donor = target_range;
            target2donor.transform(target_range_in_donor);

            // rotate search seector and find donors
            auto donor_search_sector = target_search_sector;
            target2donor.rotate(&donor_search_sector[0]);

            auto donor_list = donor_node.find_neighbors_descend(
                target_range_in_donor,
                donor_search_sector);
          
            for (auto d : donor_list) {

              auto donor_node = const_cast<block_ijk_t*>(d);
              
              // root nodes already have neighbors mapped
              if (target_node->is_root() && donor_node->is_root())
                continue;

              auto target_range_copy = target_range_in_donor;
              auto donor_range = donor_node->range();

              auto res = intersect_ranges(
                  num_dims_,
                  donor_search_sector,
                  target_range_copy,
                  donor_range);
              if (res)
              { 
                // transform back
                donor2target.transform(target_range_copy);

                target_node->add_external_neighbor(
                    *donor_node,
                    target_search_sector,
                    target_range_copy,
                    orig_neigh.ctm_me2donor_cells(),
                    orig_neigh.ctm_me2donor_verts());
              }
            } // potential donors

          } // permutations
          
        } // new targets
        
      } // sector_neighbors
    } // neighbors
  
  } // orig blocks

  return false;
   
} 

////////////////////////////////////////////////////////////////////////////////
/// Add a vertex
////////////////////////////////////////////////////////////////////////////////
int_t mesh_ijk_t::add_vertex(const real_t * vert)
{
  vertices_.insert(vertices_.end(), vert, vert+num_dims_);
  return num_vertices() - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Set a vertex
////////////////////////////////////////////////////////////////////////////////
void mesh_ijk_t::set_vertex(int_t id, const real_t * vert)
{
  std::copy(vert, vert+num_dims_, &vertices_[id*num_dims_]);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a block
////////////////////////////////////////////////////////////////////////////////
block_ijk_t & mesh_ijk_t::add_block(
    const long_t * dims,
    const int_t * verts,
    size_t nverts)
{
  // create new block
  auto bid = blocks_.size();
  blocks_.emplace_back(bid, &dims[0], &verts[0], num_dims_, nverts, *this);
  
  // add vertices->block connectivity
  std::vector<int_t> vs(&verts[0], &verts[nverts]);
  std::sort(vs.begin(), vs.end());
  block_connectivity_[std::move(vs)].emplace_back(bid);

  // add sector vertices -> block connectivity
  auto & new_block = blocks_.back();
  for (const auto & pair : new_block.verts2sector()) {
    const auto & sorted_vs = pair.first; 
    block_connectivity_[sorted_vs].emplace_back(bid);
  }


  return blocks_.back();
}

////////////////////////////////////////////////////////////////////////////////
/// Add a block
////////////////////////////////////////////////////////////////////////////////
void mesh_ijk_t::compute_neighbors()
{

  for (auto & b : blocks_) 
  {

    auto bid = b.id();
    std::set<int_t> neighbors;
  
    for (auto & pair : b.verts2sector()) {
      
      const auto & sorted_vs = pair.first;
      const auto & bs = block_connectivity_.at(sorted_vs);

      for (auto nid : bs) {
        // not me and neighbor doesnt already exist
        if (nid!=bid && !neighbors.count(nid)) {
          neighbors.emplace(nid);
          b.orient(sorted_vs, blocks_[nid]);
        }
      } // foreach connected block
    }

  }

}

////////////////////////////////////////////////////////////////////////////////
/// count leaves
////////////////////////////////////////////////////////////////////////////////
size_t mesh_ijk_t::count_leaves() const
{
  size_t n = 0;
  for (const auto & b : blocks_)
    n += b.count_leaves();
  return n;
}
  
////////////////////////////////////////////////////////////////////////////////
/// extract leaves
////////////////////////////////////////////////////////////////////////////////
void mesh_ijk_t::get_leaves(std::vector<const block_ijk_t*> & list) const
{
  for (const auto & b : blocks_)
    b.get_leaves(list);
}

std::vector<const block_ijk_t*> mesh_ijk_t::get_leaves() const
{
  std::vector<const block_ijk_t*> leaves;
  auto n = count_leaves();
  leaves.reserve(n);
  get_leaves(leaves);
  return leaves;
}

////////////////////////////////////////////////////////////////////////////////
/// Number leaves
////////////////////////////////////////////////////////////////////////////////
void mesh_ijk_t::number() {
  int_t counter = 0;
  for (auto & b : blocks_)
    b.number(counter);
}

////////////////////////////////////////////////////////////////////////////////
/// Find largest block
////////////////////////////////////////////////////////////////////////////////
const block_ijk_t * mesh_ijk_t::find_largest_block() const
{
  size_t size = 0;
  const block_ijk_t * largest = nullptr;

  for (auto & b : blocks_) {
    auto node = b.find_largest_block();
    auto nsize = node->size();
    if (nsize > size) {
      largest = node;
      size = nsize;
    }
  }

  return largest;
}

} // namespace
