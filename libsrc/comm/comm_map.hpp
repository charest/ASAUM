#ifndef COMMUNICATION_MAP_HPP
#define COMMUNICATION_MAP_HPP

#include "utils/crs.hpp"
#include "utils/errors.hpp"

#include <map>
#include <vector>

namespace prl {

struct entity_info_t;

////////////////////////////////////////////////////////////////////////////////
/// Communication map
////////////////////////////////////////////////////////////////////////////////
struct comm_map_block_t {
  int_t block_id_;

  int_t index_size_;

  std::vector<int_t> shared_users_;
  std::vector<int_t> shared_blocks_;
  crs_t<int_t, int_t> shared_ids_;
  
  std::vector<int_t> ghost_owners_;
  std::vector<int_t> ghost_blocks_;
  crs_t<int_t, int_t> ghost_ids_;

  std::map<int_t, int_t> shared_block_map_;
  std::map<int_t, int_t> ghost_block_map_;

  struct ghost_info_t {
    int_t group;
    int_t offset;
  };
  std::map<int_t, ghost_info_t> ghost_map_;

  int_t ghost_rank(int_t lid) const;
  int_t ghost_block(int_t lid) const;

  comm_map_block_t(
      int_t block_id,
      int_t index_size,
      std::vector<entity_info_t> & shared_info,
      std::vector<entity_info_t> & ghost_info);

};

////////////////////////////////////////////////////////////////////////////////
/// Communication map parent structure
////////////////////////////////////////////////////////////////////////////////
class comm_map_t {
public:
  
  using list_type = std::vector<comm_map_block_t>;
  list_type maps_;
  std::map<int_t, std::pair<int_t, list_type::iterator>> block_map_;
  
  comm_map_t(std::vector<comm_map_block_t> && maps) :
    maps_(std::move(maps))
  {
    for (size_t i=0; i<maps_.size(); ++i)
      block_map_[maps_[i].block_id_] = std::make_pair(i, maps_.begin() + i);
  }

  auto which_block(int_t i) const { return block_map_.at(i); }

  auto begin() { return maps_.begin(); }
  auto end() { return maps_.end(); }
  
  auto begin() const { return maps_.begin(); }
  auto end() const { return maps_.end(); }

  std::size_t size() const { return maps_.size(); }
  bool empty() const { return maps_.empty(); }

  auto & operator[](std::size_t i) { return maps_[i]; }
  const auto & operator[](std::size_t i) const { return maps_[i]; }

};

} // namespace

#endif
