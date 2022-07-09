#include "amr_flags.hpp"
#include "mesh_node.hpp"
#include "mesh_struct_node.hpp"

namespace prl {
  
////////////////////////////////////////////////////////////////////////////////
/// Cast as different types
////////////////////////////////////////////////////////////////////////////////
const mesh_struct_node_t * mesh_node_t::as_struct() const
{ return dynamic_cast<const mesh_struct_node_t*>(this); }

mesh_struct_node_t * mesh_node_t::as_struct()
{ return dynamic_cast<mesh_struct_node_t*>(this); }

const mesh_unstruct_node_t * mesh_node_t::as_unstruct() const
{ return dynamic_cast<const mesh_unstruct_node_t*>(this); }

mesh_unstruct_node_t * mesh_node_t::as_unstruct()
{ return dynamic_cast<mesh_unstruct_node_t*>(this); }

////////////////////////////////////////////////////////////////////////////////
/// associate a node / block together
////////////////////////////////////////////////////////////////////////////////
void mesh_node_t::set_block(std::unique_ptr<mesh_block_t> && blk)
{ 
  block_ = std::move(blk);
  block_->set_node(this);
}

void mesh_node_t::set_block(mesh_block_t * blk)
{
  block_.reset(blk);
  block_->set_node(this);
}

////////////////////////////////////////////////////////////////////////////////
/// Number leaves
////////////////////////////////////////////////////////////////////////////////
void mesh_node_t::number(int_t & counter)
{
  if (children_.size()) {
    id_ = -1;
    for (auto & c : children_)
      c ->number(counter);
  }
  else {
    id_ = counter++;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Count leaves
////////////////////////////////////////////////////////////////////////////////
size_t mesh_node_t::count_leaves() const
{
  size_t n=0;
  descend( [&](const auto &){ n++; } );
  return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Get all the leaves
////////////////////////////////////////////////////////////////////////////////
void mesh_node_t::get_leaves(std::vector<mesh_node_t*> & list) const
{
  descend([&](auto & node) {
    list.emplace_back(const_cast<mesh_node_t*>(&node));
  });
}

std::vector<mesh_node_t*> mesh_node_t::get_leaves() const
{
  auto n = count_leaves();
  std::vector<mesh_node_t*> leaves;
  leaves.reserve(n);
  get_leaves(leaves);
  return leaves;
}

////////////////////////////////////////////////////////////////////////////////
/// Get all the leaf blocks
////////////////////////////////////////////////////////////////////////////////
void mesh_node_t::get_leaf_blocks(std::vector<mesh_block_t*> & list) const
{
  descend([&](auto & node) {
    if (node.block_)
      list.emplace_back(const_cast<mesh_block_t*>(node.block_.get()));
  });
}

std::vector<mesh_block_t*> mesh_node_t::get_leaf_blocks() const
{
  auto n = count_leaves();
  std::vector<mesh_block_t*> leaves;
  leaves.reserve(n);
  get_leaf_blocks(leaves);
  return leaves;
}

////////////////////////////////////////////////////////////////////////////////
/// Free child_blocks
////////////////////////////////////////////////////////////////////////////////
void mesh_node_t::free_child_blocks()
{
  descend([&](auto & node) {
    node.block_.reset();
  });
}
  

////////////////////////////////////////////////////////////////////////////////
/// Free unused blocks
////////////////////////////////////////////////////////////////////////////////
void mesh_node_t::free_refined_blocks()
{
  if (children_.size()) {
    block_.reset();
    for (auto & c : children_)
      c->free_refined_blocks();
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Free unused blocks
////////////////////////////////////////////////////////////////////////////////
void mesh_node_t::free_marked_nodes()
{
  if (children_.size() && children_.front()->to_delete_)
    children_.clear();
  else
    for (auto & c : children_)
      c->free_marked_nodes();
}
  
////////////////////////////////////////////////////////////////////////////////
/// Part 3 - Max refinement level check
////////////////////////////////////////////////////////////////////////////////
bool mesh_node_t::level_check( amr_flags_t * flags, int_t max_level)
{
  bool is_changed = false;

  auto & bflags = flags[id_];
  auto is_rt = is_root();
  
  for (int_t dim=0; dim<num_dims_; ++dim) {
    if (level_[dim] == max_level && bflags[dim] == amr_status_t::refine) {
      bflags[dim] = amr_status_t::none;
      is_changed = true;
    } 
    else if ((is_rt || level_[dim] == 0) && bflags[dim] == amr_status_t::coarsen) {
      bflags[dim] = amr_status_t::none;
      is_changed = true;
    }
  } /// dim

  return is_changed;
}

} // namespace
