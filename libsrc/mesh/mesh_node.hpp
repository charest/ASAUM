#ifndef MESH_NODE_HPP
#define MESH_NODE_HPP

#include "mesh_block.hpp"

#include <functional>
#include <memory>
#include <vector>

namespace prl {

struct amr_flags_t;
class mesh_struct_node_t;
class mesh_unstruct_node_t;

using amr_callback_t = std::function<void(
    const amr_flags_t &,
    mesh_block_t*,
    mesh_block_t**, 
    int_t)>;

////////////////////////////////////////////////////////////////////////////////
/// The main  messh node 
////////////////////////////////////////////////////////////////////////////////
class mesh_node_t {

protected:
  
  bool to_delete_ = false;
  int_t num_dims_ = -1;
  int_t id_ = -1;
  
  std::array<int_t,3> level_ = {0, 0, 0};
  
  mesh_node_t * parent_ = nullptr;
  std::unique_ptr<mesh_block_t> block_;
  std::vector<std::unique_ptr<mesh_node_t>> children_;
  
public:

  mesh_node_t(int_t ndims, int_t id) :
    num_dims_(ndims), id_(id)
  {}

  mesh_node_t(mesh_node_t * parent) : 
    num_dims_(parent->num_dims_),
    level_{parent->level_[0], parent->level_[1],  parent->level_[2]},
    parent_(parent) 
  {}

  const mesh_struct_node_t * as_struct() const;
  mesh_struct_node_t * as_struct();
  
  const mesh_unstruct_node_t * as_unstruct() const;
  mesh_unstruct_node_t * as_unstruct();

  int_t id() const { return id_; }

  int_t num_dims() const { return num_dims_; }
  
  const int_t * level() const { return &level_[0]; }

  mesh_node_t * parent() { return parent_; }
  mesh_block_t * block() { return block_.get(); }

  mesh_node_t * child(int_t i) const { return children_[i].get(); }

  void set_unused() { id_ = -1; }
  bool is_used() const { return id_>-1; }
  bool is_root() const { return !parent_; }
  bool is_leaf() const { return children_.empty(); }

  int_t level(int_t i) const { return level_[i]; }

  bool marked_for_deletion() const { return to_delete_; }
  void mark_for_deletion() { to_delete_ = true; }
  
  const mesh_node_t * root() const { return parent_ ? parent_->root() : this; }

  bool has_children() const { return children_.size(); }
  int_t num_chidren() const { return children_.size(); }
  
  template<typename P>
  void descend(P && p) const
  {
    if (children_.size()) {
      for (auto & c : children_)
        c->descend( std::forward<P>(p) );
    }
    else {
      std::forward<P>(p)(*this);
    }
  }
  
  template<typename P>
  void descend(P && p)
  {
    if (children_.size()) {
      for (auto & c : children_)
        c->descend( std::forward<P>(p) );
    }
    else {
      std::forward<P>(p)(*this);
    }
  }

  void set_block(std::unique_ptr<mesh_block_t> && blk);
  void set_block(mesh_block_t * blk);
  
  size_t count_leaves() const;
  void number(int_t & counter);
  
  void get_leaves(std::vector<mesh_node_t*> & list) const;
  std::vector<mesh_node_t*> get_leaves() const;
  
  void get_leaf_blocks(std::vector<mesh_block_t*> & list) const;
  std::vector<mesh_block_t*> get_leaf_blocks() const;
 
  void free_refined_blocks();
  void free_child_blocks();
  void free_marked_nodes();

  virtual bool first_sibling_check(amr_flags_t *)
  { return false; }
  
  virtual bool ratio_check(amr_flags_t *, bool)
  { return false; }

  virtual bool level_check(amr_flags_t * flags, int_t max_level);
  
  virtual bool final_check_1(amr_flags_t *, const amr_flags_t *)
  { return false; }

  virtual bool final_check_2(amr_flags_t *, const amr_flags_t *)
  { return false; }
  
  virtual void change_connect(const amr_flags_t *) {}

  virtual bool refine_tree(const amr_flags_t &) { return false; }
  virtual bool coarsen_tree(const amr_flags_t &) { return false; }
  
  virtual void refine_block(const amr_flags_t &) {}
  virtual void coarsen_block() {}

  virtual void merge() {}

  virtual void clear_neighbors() {}
  virtual void find_neighbors() {}
  
  virtual ~mesh_node_t() = default;
};

////////////////////////////////////////////////////////////////////////////////
/// The main  messh node 
////////////////////////////////////////////////////////////////////////////////
class mesh_unstruct_node_t : public mesh_node_t {
public:

  mesh_unstruct_node_t(int_t ndims, int_t id=-1) : mesh_node_t(ndims, id) {}
};


} // namespace

#endif
