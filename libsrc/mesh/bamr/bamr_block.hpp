#ifndef BAMR_BLOCK_HPP
#define BAMR_BLOCK_HPP

#include "config.hpp"

#include <map>
#include <list>
#include <memory>
#include <vector>

namespace prl {

struct bamr_flags_t;

////////////////////////////////////////////////////////////////////////////////
// Binary tree
////////////////////////////////////////////////////////////////////////////////
struct bamr_block_t {
  int_t num_dims_ = -1;
  int_t id_ = -1;
  int_t label_ = -1;
  int_t split_ = -1;
  int_t sector_[3] = {0, 0, 0};
  int_t level_[3] = {0, 0, 0};
  
  bamr_block_t * parent_ = nullptr;
  bamr_block_t * left_ = nullptr;
  bamr_block_t * right_ = nullptr;

  std::map< std::array<int_t, 3>, std::vector<bamr_block_t*>> neighbors_;

  bamr_block_t(int_t id, size_t ndims)
    : num_dims_(ndims), id_(id), label_(id)
  {}
  
  bamr_block_t(int_t id, int_t label, size_t ndims)
    : num_dims_(ndims), id_(id), label_(label)
  {}

  bamr_block_t(bamr_block_t * parent) :
    sector_{parent->sector_[0], parent->sector_[1],  parent->sector_[2]},
    level_{parent->level_[0], parent->level_[1],  parent->level_[2]},
    parent_(parent)
  {
    level_[parent->split_]++;
  }

  ~bamr_block_t() {
    if (left_)  delete left_;
    if (right_) delete right_;
  }
  
  template<typename P>
  void descend(P && p) const
  {
    if (left_) {
      left_ ->descend(std::forward<P>(p));
      right_->descend(std::forward<P>(p));
    }
    else {
      std::forward<P>(p)(*this);
    }
  }

  int_t id() const { return id_; }
  int_t sector(int_t d) const { return sector_[d]; }

  size_t count_leaves() const;

  void get_leaves(std::vector<bamr_block_t*> & list) const;
  std::vector<bamr_block_t*> get_leaves() const;

  bool is_leaf() const { return !left_ || !right_; }
  bool have_children() const { return left_ && right_; }

  bool is_used() const { return id_>-1; }

  void insert_split(int_t dim);
  
  void merge()
  {
    if (have_children()) {
      left_ = nullptr;
      right_= nullptr;
      split_ = -1;
    }
  }

  bamr_block_t * sibling() const 
  { 
    if (parent_) 
      return (parent_->left_ == this) ? parent_->right_ : parent_->left_;
    else
      return nullptr;
  }


  std::vector<bamr_block_t*> & neighbors(int_t dim, int_t sector)
  {
    std::array<int_t, 3> search = {0, 0, 0};
    search[dim] = sector;
    return neighbors_.at(search);
  }


  bool same_level(bamr_block_t * other) const
  {
    return
      level_[0] == other->level_[0] &&
      level_[1] == other->level_[1] &&
      level_[2] == other->level_[2];
  }

  bamr_block_t * find_bridge(int_t dim) {
    bamr_block_t * bridge = parent_;
    while (bridge && bridge->split_ != dim)
      bridge = bridge->parent_;
    return bridge;
  }

  void lowest_level(const bamr_flags_t * flags, int_t dim, int_t & min_level);

  int_t lowest_level(const bamr_flags_t * flags, int_t dim);

  bool first_sibling_check(bamr_flags_t * flags);
  bool level_check(bamr_flags_t * flags, bool first_pass);
  bool max_level_check(bamr_flags_t * flags, int_t max_level);
  bool second_sibling_check(bamr_flags_t * flags,  const bamr_flags_t * old_flags);
  bool final_check(bamr_flags_t * flags);

  void change_connect(const bamr_flags_t * flags);

  void number(int_t & counter);
  
};


} // prl


#endif
