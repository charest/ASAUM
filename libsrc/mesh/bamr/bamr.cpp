#include "bamr.hpp"

#include "bamr_flags.hpp"
#include "math/reorder.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
// Constructor
////////////////////////////////////////////////////////////////////////////////
bamr_t::bamr_t(int_t ndims, mpi_comm_t & comm) :
  num_dims_(ndims), comm_(comm)
{
}

////////////////////////////////////////////////////////////////////////////////
// Add a root block
////////////////////////////////////////////////////////////////////////////////
bamr_block_t & bamr_t::add_block()
{
  blocks_.emplace_back(blocks_.size(), num_dims_);
  return blocks_.back();
}

bamr_block_t & bamr_t::add_block(int_t label)
{
  blocks_.emplace_back(blocks_.size(), label, num_dims_);
  return blocks_.back();
}

////////////////////////////////////////////////////////////////////////////////
/// Number leaves
////////////////////////////////////////////////////////////////////////////////
void bamr_t::number() {
  int_t counter = 0;
  for (auto & b : blocks_)
    b.number(counter);
}

////////////////////////////////////////////////////////////////////////////////
/// count leaves
////////////////////////////////////////////////////////////////////////////////
size_t bamr_t::count_leaves() const
{
  size_t n = 0;
  for (const auto & b : blocks_)
    n += b.count_leaves();
  return n;
}
  
////////////////////////////////////////////////////////////////////////////////
/// extract leaves
////////////////////////////////////////////////////////////////////////////////
void bamr_t::get_leaves(std::vector<bamr_block_t*> & list) const
{
  auto n = count_leaves();
  list.clear();
  list.reserve(n);
  for (const auto & b : blocks_)
    b.get_leaves(list);
}

std::vector<bamr_block_t*> bamr_t::get_leaves() const
{
  std::vector<bamr_block_t*> leaves;
  get_leaves(leaves);
  return leaves;
}


////////////////////////////////////////////////////////////////////////////////
// Checck
////////////////////////////////////////////////////////////////////////////////
void bamr_t::fixup(
    std::vector<bamr_block_t*> & blocks,
    std::vector<bamr_flags_t> & flags)
{
  bool first_pass = true;
  while (check(blocks, flags, first_pass))
    first_pass = false;

  change_connect(blocks, flags);
}

////////////////////////////////////////////////////////////////////////////////
// One check pass
////////////////////////////////////////////////////////////////////////////////
bool bamr_t::check(
    std::vector<bamr_block_t*> & blocks,
    std::vector<bamr_flags_t> & flags,
    bool first_pass)
{
  
  bool is_changed = false;

  while (inner_check(blocks, flags, first_pass)) {}

  if (!first_pass) {
    auto ret = final_check(blocks, flags);
    if (ret) is_changed = true;
  }

  return first_pass ? true : is_changed;
}
  
////////////////////////////////////////////////////////////////////////////////
// One check pass
////////////////////////////////////////////////////////////////////////////////
bool bamr_t::inner_check(
  std::vector<bamr_block_t*> & blocks,
  std::vector<bamr_flags_t> & flags,
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
    is_changed = is_changed || b->level_check(&flags[0], first_pass);
  
  // Part 3 - Max refinement level check
  for (auto b : blocks) 
    is_changed = is_changed || b->max_level_check(&flags[0], max_level_);

  return is_changed;
}

////////////////////////////////////////////////////////////////////////////////
// One check pass
////////////////////////////////////////////////////////////////////////////////
bool bamr_t::final_check(
    std::vector<bamr_block_t*> & blocks,
    std::vector<bamr_flags_t> & flags)
{
  bool is_changed = false;

  // store old flags
  auto old_flags = flags;
  
  // Part 4 - Second sibling check
  // 1. Levels in the direction other than the coarsening direction
  //    are the same between block and sibling.

  // check siblings in direcctinos not alligned with direction of coarsening
  for (auto b : blocks) 
    is_changed = is_changed || b->second_sibling_check(&flags[0], &old_flags[0]);

  
  // Part 5 - Impossible cases check
  for (auto b : blocks) 
    is_changed = is_changed || b->final_check(&flags[0]);

  return is_changed;
}

////////////////////////////////////////////////////////////////////////////////
// Change connectivity
////////////////////////////////////////////////////////////////////////////////
void bamr_t::change_connect(
    std::vector<bamr_block_t*> & blocks,
    const std::vector<bamr_flags_t> & flags)
{
  
  // rearrange connectivity if neihbor does not have a common parent
  for (auto b : blocks) 
    b->change_connect(&flags[0]);

}

////////////////////////////////////////////////////////////////////////////////
// Refine the mesh
////////////////////////////////////////////////////////////////////////////////
void bamr_t::refine_tree(
    std::vector<bamr_block_t*> & blocks,
    const std::vector<bamr_flags_t> & flags)
{
  for (auto & b : blocks) {
    auto bid = b->id_;
    auto & bflags = flags[bid];
    for (int_t dim=0; dim<num_dims_; ++dim)
      if (bflags[dim] == bamr_status_t::refine)
        b->insert_split(dim);
  }

}

////////////////////////////////////////////////////////////////////////////////
// Caorsen the mesh
////////////////////////////////////////////////////////////////////////////////
void bamr_t::coarsen_tree(
    std::vector<bamr_flags_t> & flags,
    int_t ndim)
{
#if 0
  for (auto & b : blocks_) {
    if (static_cast<bool>(b) && b->is_used_) {

      auto bid = b->gid_;
      auto bflags = flags[bid];
      auto parent_split = b->parent_->split_;

      while (
          bflags[0] == bamr_status_t::coarsen || 
          bflags[1] == bamr_status_t::coarsen ||
          bflags[2] == bamr_status_t::coarsen )
      {
        // coarsen along the parent split first
        if (bflags[parent_split] == bamr_status_t::coarsen) {
          b->parent_->merge(); 
          // done
          bflags[parent_split] == bamr_status_t::none;
        }
        // coarsen along other directions
        else {
          auto it = std::find(bflags.begin(), bflags.end(), bamr_status_t::coarsen);
          auto dim = std::distance(bflags.begin(), it);
          auto sector = b->sector_[dim]; 
          THROW_ERROR("Not imiplemented yet");
          // done
          bflags[dim] = bamr_status_t::none;
        }
        
      } // while

    } // used
  } // blocks
#endif
}


////////////////////////////////////////////////////////////////////////////////
// Main amr driver
////////////////////////////////////////////////////////////////////////////////
void bamr_t::amr(
    const std::vector<bamr_flags_t> & refine_flags,
    bamr_t::refine_function_t refine_mesh)
{
  
  auto blocks = get_leaves();
  auto tot_blocks = blocks.size();

  //====================================
  // Gather flags

  auto comm_rank = comm_.rank();
  auto comm_size = comm_.size();

  // gather block counts
  int_t num_blocks = refine_flags.size();
  std::vector<int_t> block_counts(comm_size);
  comm_.all_gather(num_blocks, block_counts);

  std::vector<int_t> block_offsets(comm_size+1);
  block_offsets[0] = 0;
  std::partial_sum(block_counts.begin(), block_counts.end(), &block_offsets[1]);

  // gather everyones refinement flags
  std::vector<bamr_flags_t> global_refine_flags(tot_blocks);
  comm_.all_gatherv(refine_flags, global_refine_flags, block_counts, block_offsets);

  // Make sure refinement flags wont violate any rules
  fixup(blocks, global_refine_flags);

  //====================================
  // refine blocks

  //--- refine the tree
  refine_tree(blocks, global_refine_flags);
  
  //--- number new blocks
  auto block_cnt = tot_blocks;

  // number the new blocks and track parent->child
  std::map<int_t, std::vector<bamr_block_t*>> refine_map;
  std::vector<bamr_block_t*> children;
  
  auto block_start = block_offsets[comm_rank];
  auto block_end   = block_offsets[comm_rank+1];

  for (const auto b : blocks) {
    auto bid = b->id();
    children.clear();
    b->get_leaves(children);
    // number
    for (auto c : children) {
      std::cout << "new block refined from b=" << bid <<  " sector=("
        << c->sector(0) << ", " << c->sector(1) << ", " << c->sector(2) << ")"
        << std::endl;
      c->id_ = block_cnt++;
    }
    // track mine
    if (children.size() && bid >= block_start && bid < block_end)
      refine_map[bid] = children;
  }

  //--- Create new mesh blocks
  for (auto & ref_pair : refine_map) {
    auto id = ref_pair.first;
    auto & children = ref_pair.second;
    refine_mesh(id, children.data(), children.size());
  }
  
  //====================================
  // coarsen blocks

  //coarsen_tree(global_refine_flags, ndims);

  // re-mark as used/not used
  number();
}


////////////////////////////////////////////////////////////////////////////////
// Main amr driver
////////////////////////////////////////////////////////////////////////////////
void bamr_t::refine_uniformly(int_t * levels, bamr_t::refine_function_t func)
{
  bool refine_flags[] = {false, false, false};
  bool do_refine = false;

  for (int_t d=0; d<num_dims_; ++d) {
    refine_flags[d] = levels[d] > 0;
    do_refine = do_refine || refine_flags[d];
  }

  while (do_refine) {

    refine_uniformly(refine_flags, func);

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
void bamr_t::refine_uniformly(const bool * flags, bamr_t::refine_function_t func)
{
  std::cout << "BEGIN REFINE" << std::endl;
  
  // refine in all dims
  bamr_flags_t refine_all;

  // refine in the specified dimensions
  if (flags) {
    for (int_t d=0; d<num_dims_; ++d)
      refine_all[d] = flags[d] ? bamr_status_t::refine : bamr_status_t::none;
  }
  // refine in all dimensions
  else  {
    for (int_t d=0; d<num_dims_; ++d)
      refine_all[d] = bamr_status_t::refine;
  }


  if (refine_all.has_refine()) {
    auto nblocks = count_leaves();
    std::vector<bamr_flags_t> refine_flags(nblocks, refine_all);
    amr(refine_flags, func);
  }
  
  std::cout << "DONE REFINE" << std::endl;
}

} // namespace
