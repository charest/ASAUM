
#include "entity_info.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// sort and compact entity info
////////////////////////////////////////////////////////////////////////////////
void compact(
    std::vector<entity_info_t> & entity_info,
    std::vector<int_t> & entity_ranks,
    std::vector<int_t> & entity_blocks,
    crs_t<int_t, int_t> & entity_ids)
{

  // sort the entity info
  std::sort(entity_info.begin(), entity_info.end());
  auto last = std::unique(entity_info.begin(), entity_info.end());
  entity_info.erase(last, entity_info.end());
  
  // count number of different blocks we are receiving
  int_t last_block = consts::int_max;
  int_t num_entity_blocks = 0;
  for (const auto & gi : entity_info) {
    auto blk = gi.block;
    if (blk != last_block) {
      last_block = blk;
      num_entity_blocks++;
    }
  }

  entity_ranks.clear();
  entity_blocks.clear();
  entity_ids.clear();

  entity_ranks.reserve(num_entity_blocks);
  entity_blocks.reserve(num_entity_blocks);
  entity_ids.reserve(num_entity_blocks, entity_info.size());

  auto & entity_offsets = entity_ids.offsets;
  auto & entity_indices = entity_ids.indices;

  last_block = consts::int_max;
  for (const auto & gi : entity_info) {
    auto blk = gi.block;
    if (blk != last_block) {
      entity_ranks.emplace_back(gi.rank);
      entity_blocks.emplace_back(blk);
      entity_offsets.emplace_back( entity_indices.size() );
      last_block = blk;
    }
    entity_indices.emplace_back(gi.id);
  }
    
  entity_offsets.emplace_back( entity_indices.size() );
}    

} // namespace
