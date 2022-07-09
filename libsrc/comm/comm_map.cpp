#include "comm_map.hpp"
#include "entity_info.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Communication map
////////////////////////////////////////////////////////////////////////////////
comm_map_block_t::comm_map_block_t(
    int_t block_id,
    int_t index_size,
    std::vector<entity_info_t> & shared_info,
    std::vector<entity_info_t> & ghost_info) :
  block_id_(block_id), index_size_(index_size)
{

  // sort the shared info
  compact(
      shared_info,
      shared_users_,
      shared_blocks_,
      shared_ids_);

  // pack the data
  compact(
      ghost_info,
      ghost_owners_,
      ghost_blocks_,
      ghost_ids_);

  // now build the block mappings
  int_t num_shared_blocks = shared_blocks_.size();
  for ( int_t b=0; b<num_shared_blocks; ++b) {
    auto bid = shared_blocks_[b];
    assert( shared_block_map_.count(bid) == 0 );
    shared_block_map_[bid] = b;
  }
  
  int_t num_ghost_blocks = ghost_blocks_.size();
  for ( int_t b=0; b<num_ghost_blocks; ++b) {
    auto bid = ghost_blocks_[b];
    assert( ghost_block_map_.count(bid) == 0 );
    ghost_block_map_[bid] = b;
  }
  
  // map the ghosts to the processors
  for ( int_t b=0; b<num_ghost_blocks; ++b) {
    auto start = ghost_ids_.offsets[b];
    auto end = ghost_ids_.offsets[b+1];
    for (auto i=start; i<end; ++i) {
      auto g = ghost_ids_.indices[i];
      ghost_map_[g] = {b, i};
    }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
int_t comm_map_block_t::ghost_rank(int_t lid) const
{
  auto it = ghost_map_.find(lid);
  if (it == ghost_map_.end())
    return -1;
  else
    return ghost_owners_[it->second.group];
}

////////////////////////////////////////////////////////////////////////////////
int_t comm_map_block_t::ghost_block(int_t lid) const
{
  auto it = ghost_map_.find(lid);
  if (it == ghost_map_.end())
    return -1;
  else
    return ghost_blocks_[it->second.group];
}

} // namespace
