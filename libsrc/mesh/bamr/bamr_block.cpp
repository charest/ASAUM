#include "bamr_block.hpp"

#include "bamr_flags.hpp"
#include "utils/errors.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
// Rearrange tree connectivity
////////////////////////////////////////////////////////////////////////////////
void rearrange(bamr_block_t * block1, bamr_block_t * block2, int_t split, int_t ndims)
{
  auto parent1 = block1->parent_;
  auto parent2 = block2->parent_;

  // Check levels, make sure ethey are identical.  If not, the parents
  // must be modified
  if ( !block1->same_level(block2) )
  {

    auto sibling1 = block1->sibling();
    auto sibling2 = block2->sibling();

    if (block1->split_ == -1) {
      auto parent_split = parent2->split_;
      if (block2->left_->sector_[parent_split] == -1)
        rearrange(block2->left_, sibling2->left_, parent_split, ndims);
      else
        rearrange(sibling2->left_, block2->left_, parent_split, ndims);
    }
    else if (block2->split_ == -1) {
      auto parent_split = parent1->split_;
      if (block1->left_->sector_[parent_split] == -1)
        rearrange(block1->left_, sibling1->left_, parent_split, ndims);
      else
        rearrange(sibling1->left_, block1->left_, parent_split, ndims);
    }
    else if (sibling2->split_ == -1 && sibling1->split_ != -1) {
      auto parent_split = parent1->split_;
      if (block1->left_->sector_[parent_split])
        rearrange(block1->left_, sibling1->left_, parent_split, ndims);
      else
        rearrange(sibling1->left_, block1->left_, parent_split, ndims);
    }
    else if (sibling1->split_ == -1 && sibling2->split_ != -1) {
      auto parent_split = parent2->split_;
      if (block2->left_->sector_[parent_split])
        rearrange(block2->left_, sibling2->left_, parent_split, ndims);
      else
        rearrange(sibling2->left_, block2->left_, parent_split, ndims);
    }
    else if (sibling1->split_ != -1 && sibling2->split_ != -1) {
      if (sibling1->left_->same_level(block1->left_))
      {
        auto parent_split = parent1->split_;
        if (block1->left_->sector_[parent_split] == -1)
          rearrange(block1->left_, sibling1->left_, parent_split, ndims);
        else
          rearrange(sibling1->left_, block1->left_, parent_split, ndims);
      }
      else if (sibling2->left_->same_level(block2->left_))
      {
        auto parent_split = parent2->split_;
        if (block2->left_->sector_[parent_split] == -1)
          rearrange(block2->left_, sibling2->left_, parent_split, ndims);
        else
          rearrange(sibling2->left_, block2->left_, parent_split, ndims);
      }
    }
    else {
      THROW_ERROR("No possible way to rearrange connectivity");
    }

  } // levels

  // check sector for the block pointeres parents
  if (parent1->sector_[split] != -1 && parent2->sector_[split] == 1)
    THROW_ERROR("Rearrangement cannot be performed");



  // check parents-parents
  if (!parent1->parent_ || !parent2->parent_) {
    THROW_ERROR("Already at root");
  }

  // rearrange connnectivity for parents
  if (parent1->parent_ != parent2->parent_) {
    rearrange(parent1, parent2, split, ndims);
    parent1 = block1->parent_;
    parent2 = block2->parent_;
    if (parent1->parent_ != parent2->parent_)
      THROW_ERROR("Parents still dont match");
  }

  // bridge must be found
  auto bridge = parent1->parent_;

  if (bridge->split_ != split)
    THROW_ERROR("Bridge does not have the right split");

  // change split
  auto orig_split = parent1->split_;
  bridge->split_ = orig_split;
  parent1->split_ = split;
  parent2->split_ = split;

  bamr_block_t * temp[2][2];

  // parents
  if (split != orig_split) {
    // level change
    parent1->level_[split]--;
    parent2->level_[split]--;
    parent1->level_[orig_split]++;
    parent2->level_[orig_split]++;
    // sector change
    parent1->sector_[split] = bridge->sector_[split];
    parent2->sector_[split] = bridge->sector_[split];
    parent1->sector_[orig_split] = -1;
    parent2->sector_[orig_split] = 1;
    // save children values
    auto a11 = parent1->left_ ->sector_[orig_split] == -1 ? 0 : 1;
    auto a12 = parent1->left_ ->sector_[split     ] == -1 ? 0 : 1;
    auto a21 = parent1->right_->sector_[orig_split] == -1 ? 0 : 1;
    auto a22 = parent1->right_->sector_[split     ] == -1 ? 0 : 1;

    auto b11 = parent2->left_ ->sector_[orig_split] == -1 ? 0 : 1;
    auto b12 = parent2->left_ ->sector_[split     ] == -1 ? 0 : 1;
    auto b21 = parent2->right_->sector_[orig_split] == -1 ? 0 : 1;
    auto b22 = parent2->right_->sector_[split     ] == -1 ? 0 : 1;

    temp[a11][a12] = parent1->left_;
    temp[a21][a22] = parent1->right_;
    temp[b11][b12] = parent2->left_;
    temp[b21][b22] = parent2->right_;

    // swap children pointers
    parent1->left_  = temp[a11][a12];
    parent1->right_ = temp[a21][a22];
    parent2->left_  = temp[b11][b12];
    parent2->right_ = temp[b21][b22];

    // block pointners
    parent1->left_ ->parent_ = parent1;
    parent1->right_->parent_ = parent1;
    parent2->left_ ->parent_ = parent2;
    parent2->right_->parent_ = parent2;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Number leaves
////////////////////////////////////////////////////////////////////////////////
void bamr_block_t::number(int_t & counter)
{
  if (left_) {
    id_ = -1;
    left_ ->number(counter);
    right_->number(counter);
  }
  else {
    id_ = counter++;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Count leaves
////////////////////////////////////////////////////////////////////////////////
size_t bamr_block_t::count_leaves() const
{
  size_t n=0;
  descend( [&](const auto &){ n++; } );
  return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Get all the leaves
////////////////////////////////////////////////////////////////////////////////
void bamr_block_t::get_leaves(std::vector<bamr_block_t*> & list) const
{
  auto n = count_leaves();
  list.reserve( list.size() + n );

  descend([&](auto & node) {
    list.emplace_back(const_cast<bamr_block_t*>(&node));
  });
}

std::vector<bamr_block_t*> bamr_block_t::get_leaves() const
{
  auto n = count_leaves();
  std::vector<bamr_block_t*> leaves;
  leaves.reserve(n);
  get_leaves(leaves);
  return leaves;
}
  
////////////////////////////////////////////////////////////////////////////////
// Insert a split in the tree
////////////////////////////////////////////////////////////////////////////////
void bamr_block_t::insert_split(int_t dim)
{
  if (have_children()) {
    left_ ->insert_split(dim);
    right_->insert_split(dim);
  }
  else {
    // set split dimension
    split_ = dim;
    // create new children
    left_  = new bamr_block_t(this);
    right_ = new bamr_block_t(this);
    // sectors
    left_ ->sector_[dim] = -1;
    right_->sector_[dim] =  1;
  } 
}


////////////////////////////////////////////////////////////////////////////////
// Find the lowest level
////////////////////////////////////////////////////////////////////////////////
void bamr_block_t::lowest_level(const bamr_flags_t * flags, int_t dim, int_t & min_level)
{
  if (is_used()) {
    min_level = std::min(
        min_level,
        level_[dim] + static_cast<int_t>(flags[id_][dim]));
  }
  else {
    if (left_) lowest_level(flags, dim, min_level);
    if (right_) lowest_level(flags, dim, min_level);
  }
}

int_t bamr_block_t::lowest_level(const bamr_flags_t * flags, int_t dim)
{
  auto min_level = consts::int_max;
  lowest_level(flags, dim, min_level);
  return min_level;
}

////////////////////////////////////////////////////////////////////////////////
/// Part 1 - Coarsening conflict check
///
/// 1. Sibling must have the same dimension in the direction of coarsening
/// 2. Sibling must also be flagged for coarsening.
////////////////////////////////////////////////////////////////////////////////
bool bamr_block_t::first_sibling_check( bamr_flags_t * flags)
{
  bool is_changed = false;

  //----------------------------------------------------------------------------

  auto & bflags = flags[id_];

  for (int_t dim=0; dim<num_dims_; ++dim) {
    if (bflags[dim] == bamr_status_t::coarsen) {

      auto & sector = sector_[dim];

      //--- Have siblings
      if (sector) {
        auto & neighs = neighbors(dim, sector);
        for  (auto & n : neighs) {
          auto nid = n->id_;
          auto & nflags = flags[nid];

          // neighbors not marked for coarsening
          if (n->level_[dim] != level_[dim] ||
              nflags[dim] != bamr_status_t::coarsen)
          {
            bflags[dim] = bamr_status_t::none;
            is_changed = true;
          }

          // x sibling is coarser in y, y sibling must also be flaged to coarsen in x
          else {
            for (int_t trans_dim=0; trans_dim<num_dims_; ++trans_dim) {
              if (trans_dim == dim) continue;
              if (level_[trans_dim] > n->level_[trans_dim]) {
                auto trans_sector = sector_[trans_dim];
                auto & trans_neighbors = neighbors(trans_dim, trans_sector);
                if (trans_neighbors.size() != 1 ||
                    flags[trans_neighbors[0]->id_][dim] != bamr_status_t::coarsen)
                {
                  bflags[dim] = bamr_status_t::none;
                  is_changed = true;
                }
              }
            }
          }

        } // neighbors
      
      }
      //--- No siblings, dont coarsen
      else {
        bflags[dim] = bamr_status_t::none;
        is_changed = true;
      }

    } // coarsen
  } // dim

  return is_changed;
}

////////////////////////////////////////////////////////////////////////////////
/// Part 2 - Refinement level conflict check
///
/// 1. Eliminate inadmissible cases where the level difference > 2
/// 2. Set Nothing to refine or coarsen to nothnig in some cases
////////////////////////////////////////////////////////////////////////////////
bool bamr_block_t::level_check( bamr_flags_t * flags, bool first_pass)
{
  bool is_changed = false;
  
    
  auto & bflags = flags[id_];
  
  for (int_t dim=0; dim<num_dims_; ++dim) {
    if (bflags[dim] != bamr_status_t::none) {

      bool to_change = false;
  
      //--- Loop through all neighbors
      for (auto & sector_neighbors : neighbors_) {
        for (auto & n : sector_neighbors.second) {

          auto nid = n->id_;
          auto & nflags = flags[nid];
    
          auto leveldiff =
            (level_[dim] + static_cast<int_t>(bflags[dim])) -
            (n->level_[dim] + static_cast<int_t>(nflags[dim]));

          // if the coarser neighbor flagged coarsen, and the finer neighbor
          // flagged refine, 
          if (std::abs(leveldiff) > 2) {
            nflags[dim] = bamr_status_t::none;
            to_change = true;
          }
          else if (!first_pass) {
            if (leveldiff < -1) {
              nflags[dim] = (nflags[dim] == bamr_status_t::none) ?
                bamr_status_t::refine : bamr_status_t::none;
              is_changed = true;
            }
            else if (leveldiff > 1) {
              to_change = true;
            }
          }

        } // neighbors

      } // sector

      if (to_change) {
        bflags[dim] = bamr_status_t::none;
        is_changed = true;
      }

    } // refine
  } // dim

  return is_changed;
 
}

////////////////////////////////////////////////////////////////////////////////
/// Part 3 - Max refinement level check
////////////////////////////////////////////////////////////////////////////////
bool bamr_block_t::max_level_check( bamr_flags_t * flags, int_t max_level)
{
  bool is_changed = false;

  auto & bflags = flags[id_];
  
  for (int_t dim=0; dim<num_dims_; ++dim) {
    if (level_[dim] == max_level && bflags[dim] == bamr_status_t::refine) {
      bflags[dim] = bamr_status_t::none;
      is_changed = true;
    } 
  } /// dim

  return is_changed;
}

////////////////////////////////////////////////////////////////////////////////
/// Part 4 - Second sibling check
///
/// 1. Levels in the direction other than the coarsening direction
///    are the same between block and sibling.
////////////////////////////////////////////////////////////////////////////////
bool bamr_block_t::second_sibling_check(
    bamr_flags_t * flags,
    const bamr_flags_t * old_flags)
{
  bool is_changed = false;

  auto bid = id_;
  auto & bflags = flags[bid];

  for (int_t dim=0; dim<num_dims_; ++dim) {
    if (bflags[dim] == bamr_status_t::coarsen) {

      auto & sector = sector_[dim];
      auto & neighs = neighbors(dim, sector);
      
      for (auto n : neighs)
      {
        auto nid = n->id_;
        auto & nflags = flags[nid];
          
        //--- two dimensions
        for (int_t trans_dim=0; trans_dim<num_dims_ && num_dims_==2; ++trans_dim) {
          if (trans_dim == dim) continue;
        
          auto transflags = old_flags[bid][trans_dim];
          auto neighflags = old_flags[nid][trans_dim];
          auto leveldiff = 
            (n->level_[trans_dim] + static_cast<int_t>(neighflags)) -
            (   level_[trans_dim] + static_cast<int_t>(transflags));

          if (leveldiff != 0) {
            // block set too coarse, but can be forced to refine
            if (leveldiff == 1 && transflags == bamr_status_t::none) {
              bflags[trans_dim] = bamr_status_t::refine;
            }
            // neighbor block is too coarsee, but can be forced to refine
            else if (leveldiff == -1 && neighflags == bamr_status_t::none) {
              nflags[trans_dim] = bamr_status_t::refine;
            }
            else {
              bflags[dim] = bamr_status_t::none;
              nflags[dim] = bamr_status_t::none;
            }

            is_changed = true;
          } // leveldiff

        } // trans_dim
        
        //--- three dimensions
        for (int_t trans_dim=0; trans_dim<num_dims_ && num_dims_==3; ++trans_dim) {
          if (trans_dim == dim) continue;

          auto other_dim = num_dims_ - dim - trans_dim;
        
          auto transflags1 = old_flags[bid][trans_dim];
          auto neighflags1 = old_flags[nid][trans_dim];
          auto leveldiff1 = 
            (n->level_[trans_dim] + static_cast<int_t>(neighflags1)) -
            (   level_[trans_dim] + static_cast<int_t>(transflags1));
          
          auto transflags2 = old_flags[bid][other_dim];
          auto neighflags2 = old_flags[nid][other_dim];
          auto leveldiff2 = 
            (n->level_[other_dim] + static_cast<int_t>(neighflags2)) -
            (   level_[other_dim] + static_cast<int_t>(transflags2));

          if (leveldiff1 != 0) {
            // block set too coarse, but can be forced to refine
            if (leveldiff1 == 1 && transflags1 == bamr_status_t::none) {
              // block is too coarse in y but can be forced to refine
              // (since no level diff in z)
              if (leveldiff2 == 0) {
                bflags[trans_dim] = bamr_status_t::refine;
              }
              // block is too coarse in y and z, but can be forced to refine in each
              else if (leveldiff2 == 1 && transflags2 == bamr_status_t::none) {
                bflags[trans_dim] = bamr_status_t::refine;
                bflags[other_dim] = bamr_status_t::refine;
              }
              // nothing can be done, abort mission
              else {
                bflags[dim] = bamr_status_t::none;
                nflags[dim] = bamr_status_t::none;
              }
            }
            // neighbor block is too coarse in y but can be forced to refine
            else if (leveldiff1 == -1 && neighflags1 == bamr_status_t::none) {
              // block is too coarse in y but can be forced to refine
              // (since no level diff in z)
              if (leveldiff2 == 0) {
                neighflags1 = bamr_status_t::refine;
              }
              // block is too coarse in y and z, but can be forced to refine in each
              else if( leveldiff2 == 1 && transflags2 == bamr_status_t::none) {
                bflags[trans_dim] = bamr_status_t::refine;
                bflags[other_dim] = bamr_status_t::refine;
              }
              // nothing can be done, abort mission
              else {
                bflags[dim] = bamr_status_t::none;
                nflags[dim] = bamr_status_t::none;
              }
            }

            is_changed = true;
          
          } // leveldiff

        } // trans_dim

      } // neighbors

    } // coarsen
  } // dim

  return is_changed;
}

////////////////////////////////////////////////////////////////////////////////
/// Part 5 - Impossible cases check
////////////////////////////////////////////////////////////////////////////////
bool bamr_block_t::final_check( bamr_flags_t * flags)
{
  bool is_changed = false;
      
  auto & bflags = flags[id_];
  auto & blevel = level_;

  for (int_t dim=0; dim<num_dims_; ++dim) {
    auto & bflag = bflags[dim];
    if (bflag == bamr_status_t::coarsen) {

      auto parent = parent_;
      bamr_block_t * bridge = find_bridge(dim);

      if (!bridge) {
        bflag = bamr_status_t::none;
        is_changed = true;
      }
      else if (parent != bridge) {
        // if the lowest level of refinement in the other dimensions
        // are lower than this block, do not allow coarsening
        bool change = true;
        for (int_t trans_dim=0; trans_dim<num_dims_; ++trans_dim) {
          if (trans_dim == dim) continue;
          auto lowest_level = bridge->lowest_level(&flags[0], trans_dim);
          if  (lowest_level >= blevel[trans_dim] + static_cast<int_t>(bflags[trans_dim]))
            change = false;
        }
        if (change) {
          bflag = bamr_status_t::none;
          is_changed = true;
        }
      }

    }
  } // dim

  return is_changed;
}
      
////////////////////////////////////////////////////////////////////////////////
/// Change connectivity
////////////////////////////////////////////////////////////////////////////////
void bamr_block_t::change_connect(const bamr_flags_t * flags)
{
  auto & bflags = flags[id_];

  for (int_t dim=0; dim<num_dims_; ++dim) {
    if (bflags[dim] == bamr_status_t::coarsen) {
    
      auto & sector = sector_[dim];
      auto & neighs = neighbors(dim, sector);

      //--- Only one neighbor
      if (neighs.size()==1) {

        // if neighbor has same resolution in all dims
        auto & neigh = neighs[0];
        if (same_level(neigh))
        {
          for (int_t trans_dim=0; trans_dim>num_dims_; ++trans_dim) {
            if (trans_dim==dim) continue;
            // if block is not also flagged to coarsen in transverse direction,
            // make sure parent is split in normal direction
            bool option1 = 
              (bflags[trans_dim] != bamr_status_t::coarsen) &&
              (parent_->split_ != dim);
            // if block is flagged to coarsen in tranverse direction,
            // but both parent and grandparent are split in transversee direction
            bool option2 =
              (bflags[trans_dim] == bamr_status_t::coarsen) && 
              (parent_->split_ == trans_dim) &&
              (parent_->parent_->split_ == trans_dim);
            // if either true, rearrange to place a split
            if (option1 || option2) {
              if (sector == -1)
                rearrange(this, neigh, dim, num_dims_);
              else
                rearrange(neigh, this, dim, num_dims_);
            } // options

          } // trans_dim
        } // same refinement level

      } 
      //--- Two neighbors
      else if (neighs.size()==2) {
        if (parent_->split_ != dim) {
          auto & neigh = neighs[0];
          if (sector == -1)
            rearrange(this, neigh->parent_, dim, num_dims_);
          else
            rearrange(neigh->parent_, this, dim, num_dims_);
        }
      }

    } // coarseen
  } // dim

}

} // namespace
