#ifndef ENTITY_INFO_HPP
#define ENTITY_INFO_HPP

#include "config.hpp"

#include "utils/crs.hpp"

#include <vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// A struct for entity info
////////////////////////////////////////////////////////////////////////////////
struct entity_info_t {
  int_t id;
  int_t rank;
  int_t block;
  int_t sort;
  
  entity_info_t(int_t i, int_t r, int_t blk) :
    id(i), rank(r), block(blk), sort(i) {}

  entity_info_t(int_t i, int_t r, int_t blk, int_t srt) :
    id(i), rank(r), block(blk), sort(srt) {}

  bool operator<(const entity_info_t & b) {
    if (rank != b.rank) return (rank < b.rank);
    if (block != b.block) return (block < b.block);
    if (sort != b.sort) return  (sort < b.sort);
    return false;
  }
  
  bool operator==(const entity_info_t & b) {
    if (rank != b.rank) return false;
    if (block != b.block) return false;
    return (sort == b.sort);
  }
};

/// sort and compact entity info
void compact(
    std::vector<entity_info_t> & entity_info,
    std::vector<int_t> & entity_ranks,
    std::vector<int_t> & entity_blocks,
    crs_t<int_t, int_t> & entity_ids);

} // namespace

#endif
