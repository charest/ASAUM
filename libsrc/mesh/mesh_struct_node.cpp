#include "mesh_struct_node.hpp"
#include "mesh_struct_block.hpp"

#include "amr_flags.hpp"

#include <cmath>
#include <iomanip>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// scale a range
////////////////////////////////////////////////////////////////////////////////
void scale_range(
    int_t ndims,
    const int_t * src_level,
    const int_t * tgt_level,
    range_ijk_t & range)
{
  for (int_t d=0; d<ndims; ++d) {
    auto delta = tgt_level[d] - src_level[d];
    if (delta < 0)
      range.div(d, std::pow(2, -delta));
    else if (delta > 0)
      range.mult(d, std::pow(2, delta));
  }
}

////////////////////////////////////////////////////////////////////////////////
/// scale a range
////////////////////////////////////////////////////////////////////////////////
bool do_ranges_intersect(
    const int_t nd,
    const sector_t & sector,
    const int_t * my_level,
    const range_ijk_t & my_range,
    const int_t * neigh_level,
    const range_ijk_t & neigh_range)
{
  for (int_t dim=0; dim<nd; ++dim) {

    int_t my_mult = 1;
    int_t nb_mult = 1;
    if (my_level[dim] != neigh_level[dim]) {
      auto delta = my_level[dim] - neigh_level[dim];
      if (delta<0) {
        my_mult = std::pow(2, -delta);
      }
      else {
        nb_mult = std::pow(2, delta);
      }
    }

    auto my_range_begin = my_mult * my_range.begin[dim];
    auto my_range_end   = my_mult * my_range.end[dim];
    
    auto neigh_range_begin = nb_mult * neigh_range.begin[dim];
    auto neigh_range_end   = nb_mult * neigh_range.end[dim];

    if (sector[dim] == 0) {
      auto beg = std::max(my_range_begin, neigh_range_begin);
      auto end = std::min(my_range_end, neigh_range_end);
      auto len = end - beg;
      if (len<1) return false;
    }
    else {
      if (sector[dim] < 0) {
        if (neigh_range_begin >= my_range_begin ||
            neigh_range_end   <  my_range_begin)
          return false;
      }
      else /*if (sector[dim] > 0)*/ {
        if (neigh_range_begin >  my_range_end ||
            neigh_range_end   <= my_range_end)
          return false;
      }
    } // else
  } // dim

  return true;
}


bool do_ranges_touch_left(
    int_t dim,
    const int_t * level_a,
    const int_t * level_b,
    const range_ijk_t & range_a,
    const range_ijk_t & range_b)
{
  int_t mult_a = 1;
  int_t mult_b = 1;
  if (level_a[dim] != level_b[dim]) {
    auto delta = level_a[dim] - level_b[dim];
    if (delta<0) {
      mult_a = std::pow(2, -delta);
    }
    else {
      mult_b = std::pow(2, delta);
    }
  }
  return (
      mult_b*range_b.begin[dim] <  mult_a*range_a.begin[dim] &&
      mult_b*range_b.end[dim]   >= mult_a*range_a.begin[dim] );
}

bool do_ranges_touch_right(
    int_t dim,
    const int_t * level_a,
    const int_t * level_b,
    const range_ijk_t & range_a,
    const range_ijk_t & range_b)
{
  int_t mult_a = 1;
  int_t mult_b = 1;
  if (level_a[dim] != level_b[dim]) {
    auto delta = level_a[dim] - level_b[dim];
    if (delta<0) {
      mult_a = std::pow(2, -delta);
    }
    else {
      mult_b = std::pow(2, delta);
    }
  }
  return (
      mult_b*range_b.begin[dim] <= mult_a*range_a.end[dim] &&
      mult_b*range_b.end[dim]   >  mult_a*range_a.end[dim] );
}

bool intersect_ranges(
    const int_t nd,
    const sector_t & sector,
    const int_t * my_level,
    range_ijk_t & my_range,
    const int_t * neigh_level,
    range_ijk_t & neigh_range)
{
  for (int_t dim=0; dim<nd; ++dim) {
    
    int_t my_mult = 1;
    int_t nb_mult = 1;
    if (my_level[dim] != neigh_level[dim]) {
      auto delta = my_level[dim] - neigh_level[dim];
      if (delta<0) {
        my_mult = std::pow(2, -delta);
      }
      else {
        nb_mult = std::pow(2, delta);
      }
    }

    std::cout << "mults " << my_mult << " " << nb_mult << std::endl;
    
    auto my_range_begin = my_mult * my_range.begin[dim];
    auto my_range_end   = my_mult * my_range.end[dim];
    
    auto neigh_range_begin = nb_mult * neigh_range.begin[dim];
    auto neigh_range_end   = nb_mult * neigh_range.end[dim];

    if (sector[dim] == 0) {
      auto beg = std::max(my_range_begin, neigh_range_begin);
      auto end = std::min(my_range_end, neigh_range_end);
      auto len = end - beg;
      neigh_range.begin[dim] = beg / nb_mult;
      my_range   .begin[dim] = beg / my_mult;
      neigh_range.end  [dim] = end / nb_mult;
      my_range   .end  [dim] = end / my_mult;
      if (len<1) return false;
    }
    else {
      if (sector[dim] < 0) {
        auto end = std::min( my_range_begin, neigh_range_end );
        my_range.begin[dim] = end / my_mult;
        neigh_range.end[dim] = end / nb_mult;
      }
      else if (sector[dim] > 0) {
        auto beg = std::max( my_range_end, neigh_range_begin);
        neigh_range.begin[dim] = beg / nb_mult;
        my_range.end[dim] =  beg / my_mult;
      }
      auto len = neigh_range.end[dim] - neigh_range.begin[dim];
      if (len<1) return false;
    }
  }

  return true;
}


////////////////////////////////////////////////////////////////////////////////
/// constructors
////////////////////////////////////////////////////////////////////////////////
mesh_struct_node_t::mesh_struct_node_t(
    const long_t * dims,
    int_t ndims,
    int_t id) :
  mesh_node_t(ndims, id),
  range_(num_dims_, dims)
{ std::copy_n(dims, ndims, &dims_[0]); }

mesh_struct_node_t::mesh_struct_node_t(
    mesh_struct_node_t * parent,
    int_t split,
    int_t sector,
    long_t midp) :
  mesh_node_t(parent),
  dims_{parent->dims_[0], parent->dims_[1], parent->dims_[2]},
  range_(parent_->as_struct()->range_)
{
  level_[split]++;
  range_.begin[split] *= 2;
  range_.end[split] *= 2;
  if (sector < 0)
    range_.end[split] = midp*2;
  else if (sector > 0)
    range_.begin[split] = midp*2;
  dims_[split] = range_.size(split);
}

////////////////////////////////////////////////////////////////////////////////
// Insert a split in the tree
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::insert_split(int_t dim)
{
  if (children_.size()) {
    children_[0]->as_struct()->insert_split(dim);
    children_[1]->as_struct()->insert_split(dim);
  }
  else {
    // set split dimension
    split_ = dim;
    // create new children
    auto midp = range_.begin[dim] + dims_[dim] / 2;
    auto left  = std::make_unique<mesh_struct_node_t>(this, dim, -1, midp);
    auto right = std::make_unique<mesh_struct_node_t>(this, dim,  1, midp);
    // add them to our list
    children_.reserve(2);
    children_.emplace_back(std::move(left));
    children_.emplace_back(std::move(right));
  } 
}

////////////////////////////////////////////////////////////////////////////////
// Merge the tree
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::merge()
{
  std::cout << "markinng " << children_.size() << " for delletion" << std::endl;
  for (auto & c : children_) {
    std::cout << "marked " << c->id() << " for deletion" << std::endl;
    c->mark_for_deletion();
  }
  split_ = -1;
}
  
////////////////////////////////////////////////////////////////////////////////
/// Part 1 - Coarsening conflict check
///
/// 1. Sibling must have the same dimension in the direction of coarsening
/// 2. Sibling must also be flagged for coarsening.
////////////////////////////////////////////////////////////////////////////////
bool mesh_struct_node_t::first_sibling_check(amr_flags_t * flags)
{
  bool is_changed = false;

  auto & bflags = flags[id_];
  std::vector<const mesh_struct_node_t *> siblings;

  for (int_t dim=0; dim<num_dims_; ++dim) {
    if (bflags[dim] == amr_status_t::coarsen) {

      siblings.clear();
      find_siblings(dim, siblings);

      //--- Have siblings
      if (siblings.size()) {
        
        for  (auto & n : siblings) {
          auto nid = n->id_;
          auto & nflags = flags[nid];

          // neighbors not marked for coarsening
          if (n->level_[dim] != level_[dim] ||
              nflags[dim] != amr_status_t::coarsen)
          {
            bflags[dim] = amr_status_t::none;
            is_changed = true;
          }

        } // siblings
      
      }
      //--- No siblings, dont coarsen
      else {
        bflags[dim] = amr_status_t::none;
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
bool mesh_struct_node_t::ratio_check(
    amr_flags_t * flags,
    bool first_pass)
{
  bool is_changed = false;
  
  auto & bflags = flags[id_];
  
  for (int_t dim=0; dim<num_dims_; ++dim) {
    if (bflags[dim] != amr_status_t::none) {

      //--- Loop through all neighbors
      for (auto & sector_neighbors : neighbors_) {
        for (auto & n : sector_neighbors.second) {

          auto nid = n.id();
          auto & nflags = flags[nid];
          const auto & nlevel = n.node()->level();
    
          auto leveldiff =
            (nlevel[dim] + static_cast<int_t>(nflags[dim])) -
            (level_[dim] + static_cast<int_t>(bflags[dim]));
          std::cout << "me id " << id_
            << " " << std::setw(2) <<  level_[dim] << " " << std::setw(2) <<
            static_cast<int_t>(bflags[dim]) << 
            " neigh id " << nid << " "  << std::setw(2) << nlevel[dim] << " " << std::setw(2) << 
            static_cast<int_t>(nflags[dim]) <<
            " leveldiff " << leveldiff << " now ";
          
          // if the coarser neighbor flagged coarsen, and the finer neighbor
          // flagged refine, 
          if (std::abs(leveldiff) > 2) {
            bflags[dim] = amr_status_t::none;
            nflags[dim] = amr_status_t::none;
            std::cout << "both changed ";
            is_changed = true;
          }
          else if (!first_pass) {
            if (leveldiff < -1) {
              std::cout << "neigh changed ";
              nflags[dim] = (nflags[dim] == amr_status_t::none) ?
                amr_status_t::refine : amr_status_t::none;
              is_changed = true;
            }
            else if (leveldiff > 1) {
              std::cout << "me changed ";
              bflags[dim] = (bflags[dim] == amr_status_t::none) ?
                amr_status_t::refine : amr_status_t::none;
              is_changed = true;
            }
          }
          
          leveldiff =
            (nlevel[dim] + static_cast<int_t>(nflags[dim])) -
            (level_[dim] + static_cast<int_t>(bflags[dim]));
          std::cout << "leveldiff " << leveldiff << std::endl;

        } // neighbors

      } // sector

    } // refine
  } // dim

  return is_changed;
 
}


////////////////////////////////////////////////////////////////////////////////
/// Return a structured block
////////////////////////////////////////////////////////////////////////////////
mesh_struct_block_t * mesh_struct_node_t::struct_block()
{ return dynamic_cast<mesh_struct_block_t*>(block()); }

////////////////////////////////////////////////////////////////////////////////
/// Determine block range
////////////////////////////////////////////////////////////////////////////////
range_ijk_t mesh_struct_node_t::range(const sector_t & sector) const
{
  range_ijk_t r;

  auto ndims = num_dims();
  for (int_t d=0; d<ndims; ++d) {
    if (sector[d] == 0) {
      r.begin[d] = 0;
      r.end[d] = dims_[d];
    }
    else if (sector[d] < 0) {
      r.begin[d] = 0;
      r.end[d] = 1;
    }
    else {
      r.begin[d] = dims_[d]-1;
      r.end[d] = dims_[d];
    }
  }

  return r;
}

range_ijk_t mesh_struct_node_t::range_wrt_root(const sector_t & sector) const
{
  range_ijk_t r;

  auto ndims = num_dims();
  for (int_t d=0; d<ndims; ++d) {
    if (sector[d] == 0) {
      r.begin[d] = range_.begin[d];
      r.end[d] = range_.end[d];
    }
    else if (sector[d] < 0) {
      r.begin[d] = range_.begin[d];
      r.end[d] = range_.begin[d]+1;
    }
    else {
      r.begin[d] = range_.end[d]-1;
      r.end[d] = range_.end[d];
    }
  }

  return r;
}


////////////////////////////////////////////////////////////////////////////////
/// Make a range relative to me
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::make_relative(long_t * pos) const
{
  for (int_t dim=0; dim<num_dims_; ++dim)
    pos[dim] -= range_.begin[dim];
}

void mesh_struct_node_t::make_relative(range_ijk_t & rng) const
{
  for (int_t dim=0; dim<num_dims_; ++dim) {
    rng.begin[dim] -= range_.begin[dim];
    rng.end[dim]   -= range_.begin[dim];
  }
}

////////////////////////////////////////////////////////////////////////////////
/// clear neighbor info
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::clear_neighbors()
{
  neighbors_.clear();
  shared_.clear();

  for (auto & c : children_)
    c->clear_neighbors();
}

////////////////////////////////////////////////////////////////////////////////
/// add an internal neighbor
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::add_internal_neighbor(
    mesh_struct_node_t * neigh,
    const sector_t & my_sector)
{
  std::cout << "adding internal neighbor between " << id_ << " and " << neigh->id() << std::endl;
  std::cout << "my level " << level_[0] << " vs " << neigh->level_[0] << std::endl;
  assert(std::abs(neigh->level_[0] - level_[0]) < 2);

  //--- Get ranges
  auto my_range = range_wrt_root(my_sector);
  auto donor_range = neigh->range();
  
  scale_range(num_dims_, &neigh->level_[0], &level_[0], donor_range);

  std::cout << "my range " << my_range << std::endl;
  std::cout << "donor range " << donor_range << std::endl;
  //--- Intersect with neighbor
  auto res = intersect_ranges(
      num_dims_,
      my_sector,
      my_range,
      donor_range);
  assert(res);
  UNUSED(res);
  
  std::cout << "res " << res << std::endl;
  std::cout << "int my range " << my_range << std::endl;
  std::cout << "int donor range " << donor_range << std::endl;
  
  auto scaled_range = my_range;
  scale_range(num_dims_, &level_[0], &neigh->level_[0], scaled_range);
  
  //--- Compute a known position in target frame
  std::array<long_t,3> my_cell_pos{0, 0, 0}, my_vert_pos{0, 0, 0};
  std::array<long_t,3> neigh_cell_pos{0, 0, 0}, neigh_vert_pos{0, 0, 0};
  for (int_t d=0; d<num_dims_; ++d) {
    if (my_sector[d] < 0)
      scaled_range.end[d] = scaled_range.begin[d]+1;
    else if (my_sector[d] > 0)
      scaled_range.begin[d] = scaled_range.end[d]-1;
    my_cell_pos[d] = my_range.begin[d] + my_sector[d];
    my_vert_pos[d] = my_sector[d]>0 ? my_range.begin[d]+1 : my_range.begin[d];
    neigh_cell_pos[d] = scaled_range.begin[d] + my_sector[d];
    neigh_vert_pos[d] = my_sector[d]>0 ? scaled_range.begin[d]+1 : scaled_range.begin[d];
  }
  
  //--- make positions relative
  make_relative(&my_cell_pos[0]);
  make_relative(&my_vert_pos[0]);

  neigh->make_relative(&neigh_cell_pos[0]);
  neigh->make_relative(&neigh_vert_pos[0]);
  
  //--- compute transformation
  ctm_t me2neigh_cells(num_dims_), me2neigh_verts(num_dims_);
  
  //--- compute offsets
  for (int_t d=0; d<num_dims_; ++d) {
    me2neigh_cells.offset(d) = neigh_cell_pos[d] - my_cell_pos[d];
    me2neigh_verts.offset(d) = neigh_vert_pos[d] - my_vert_pos[d];
  }
  
  //--- Make ranges relative
  make_relative(my_range);
  
  //--- Add my neighbor 
  add_neighbor(my_sector, neigh, me2neigh_cells, me2neigh_verts, my_range);
  std::cout << "neighbor range " << my_sector << " " << my_range << std::endl;

  //--- Add shared info
  for (int_t d=0; d<num_dims_; ++d) {
    my_range.begin[d] += my_sector[d];
    my_range.end[d] += my_sector[d];
  }
  me2neigh_cells.transform(my_range);

  neigh->add_shared(this, my_sector, my_range);
  std::cout << "neighbor range " << my_sector << " " << my_range << std::endl;

  std::cout << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a block neighbor
/// \note neigh_range is in the donor frame
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::add_external_neighbor(
    mesh_struct_node_t * neigh,
    sector_t my_sector,
    transform_t me2neigh_cells,
    transform_t me2neigh_verts)
{
  std::cout << "adding externa neighbor between " << id_ << " and " << neigh->id() << std::endl;
  std::cout << "my level " << level_[0] << " vs " << neigh->level_[0] << std::endl;
  assert(std::abs(neigh->level_[0] - level_[0]) < 2);
  
  //--- Find donnor sector
  auto donor_sector = flip(my_sector);
  me2neigh_cells.rotate(&donor_sector[0]);

  //--- Get ranges
  auto my_range = range_wrt_root(my_sector);
  auto donor_range = neigh->range();
  
  std::cout << "my range " << my_range << std::endl;
  std::cout << "donor range " << donor_range << std::endl;

  //--- intersect
  auto neigh2me = reverse(me2neigh_cells);
  neigh2me.scale1(&level_[0]);
  scale_range(num_dims_, &neigh->level_[0], &level_[0], donor_range);
  neigh2me.transform(donor_range);
  
  std::cout << "transformed donor range " << donor_range << std::endl;

  auto res = intersect_ranges(
      num_dims_,
      my_sector,
      my_range,
      donor_range);
  assert(res);
  UNUSED(res);

  std::cout << "res " << res << std::endl;
  std::cout << "int my range " << my_range << std::endl;
  std::cout << "int donor range " << donor_range << std::endl;
  
  auto scaled_range = my_range;
  scale_range(num_dims_, &level_[0], &neigh->level_[0], scaled_range);

  //--- Compute a known position in target frame
  std::array<long_t,3> my_cell_pos{0, 0, 0}, my_vert_pos{0, 0, 0};
  std::array<long_t,3> neigh_cell_pos{0, 0, 0}, neigh_vert_pos{0, 0, 0};
  for (int_t d=0; d<num_dims_; ++d) {
    if (my_sector[d] < 0)
      scaled_range.end[d] = scaled_range.begin[d]+1;
    else if (my_sector[d] > 0)
      scaled_range.begin[d] = scaled_range.end[d]-1;
    my_cell_pos[d] = my_range.begin[d] + my_sector[d];
    my_vert_pos[d] = my_sector[d]>0 ? my_range.begin[d]+1 : my_range.begin[d];
    neigh_cell_pos[d] = scaled_range.begin[d] + my_sector[d];
    neigh_vert_pos[d] = my_sector[d]>0 ? scaled_range.begin[d]+1 : scaled_range.begin[d];
  }

  //--- transform it to donor frame
  std::cout << "my_cell pos (" << my_cell_pos[0] << ") vs "
    << " neig pos  (" << neigh_cell_pos[0] << ")" << std::endl;
  std::cout << "my_vert pos (" << my_vert_pos[0] << ") vs "
    << " neig pos  (" << neigh_vert_pos[0] << ")" << std::endl;

  me2neigh_cells.scale1(&neigh->level_[0]);
  me2neigh_verts.scale2(&neigh->level_[0]);
  me2neigh_cells.transform(&neigh_cell_pos[0]);
  me2neigh_verts.transform(&neigh_vert_pos[0]);

  std::cout << "my_cell pos (" << my_cell_pos[0] << ") vs "
    << " neig pos  (" << neigh_cell_pos[0] << ")" << std::endl;
  std::cout << "my_vert pos (" << my_vert_pos[0] << ") vs "
    << " neig pos  (" << neigh_vert_pos[0] << ")" << std::endl;

  //--- make positions relative
  make_relative(&my_cell_pos[0]);
  make_relative(&my_vert_pos[0]);

  neigh->make_relative(&neigh_cell_pos[0]);
  neigh->make_relative(&neigh_vert_pos[0]);
  
  //--- Make initial transformation
  me2neigh_cells.zero_offsets();
  me2neigh_verts.zero_offsets();

  auto new_neigh_cell_pos = my_cell_pos;
  auto new_neigh_vert_pos = my_vert_pos;

  me2neigh_cells.transform(&new_neigh_cell_pos[0]);
  me2neigh_verts.transform(&new_neigh_vert_pos[0]);

  //--- compute offsets
  for (int_t d=0; d<num_dims_; ++d) {
    me2neigh_cells.offset(d) = neigh_cell_pos[d] - new_neigh_cell_pos[d];
    me2neigh_verts.offset(d) = neigh_vert_pos[d] - new_neigh_vert_pos[d];
  }

  //--- Make ranges relative
  make_relative(my_range);
  
  //--- Add targets neighbor
  add_neighbor(
      my_sector,
      neigh,
      compact(me2neigh_cells),
      compact(me2neigh_verts),
      my_range);
  std::cout << "neighbor range " << my_sector << " " << my_range << std::endl;
  
  //--- Add shared info
  for (int_t d=0; d<num_dims_; ++d) {
    my_range.begin[d] += my_sector[d];
    my_range.end[d] += my_sector[d];
  }
  me2neigh_cells.rotate(&my_sector[0]);
  me2neigh_cells.transform(my_range);
  neigh->add_shared(this, my_sector, my_range);
  
  std::cout << "neighbor range " << my_sector << " " << my_range << std::endl;

  std::cout << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Find my sibling
////////////////////////////////////////////////////////////////////////////////
mesh_struct_node_t * mesh_struct_node_t::sibling() const
{
  if (parent_)
    return
      (parent_->child(0) == static_cast<const mesh_struct_node_t*>(this)) ?
        parent_->child(1)->as_struct() : parent_->child(0)->as_struct();
  else
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Find left and right nodes
////////////////////////////////////////////////////////////////////////////////
mesh_struct_node_t * mesh_struct_node_t::left(const mesh_struct_node_t * c) const
{
  if (c == child(1))
   return child(0)->as_struct();
  else
    return nullptr;
}

mesh_struct_node_t * mesh_struct_node_t::right(const mesh_struct_node_t * c) const
{
  if (c == child(0))
   return child(1)->as_struct();
  else
    return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Find my sibling
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::find_siblings(
    int_t dim,
    std::vector<const mesh_struct_node_t *> & siblings) const
{
  std::cout << "OOKING FOR SIBLINGS " << id() << std::endl;
  if (parent_) {

    auto parent = parent_->as_struct();

    // split is in the same directio
    if (parent->split_ == dim) {

      auto l = parent->child(0)->as_struct();
      auto r = parent->child(1)->as_struct();
      std::cout << "l " << l << " r " << r << " " << this << std::endl;
      // go down right side, stay left
      if (l == this) {
        std::cout << "going down right" << std::endl;
        r->find_siblings_descend(dim, 0, siblings);
      }
      // go down left side, stay right
      else {
        std::cout << "going down left" << std::endl;
        l->find_siblings_descend(dim, 1, siblings);
      }

    }
    // split is not in same direction, keep goiing
    else {
      parent->find_siblings(dim, siblings);
    }
  }
  std::cout << "DONE" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Find my sibling
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::find_siblings_descend(
    int_t dim,
    bool go_right,
    std::vector<const mesh_struct_node_t *> & siblings) const
{
  if (children_.empty()) {
    std::cout << "addingn " << id() << std::endl;
    siblings.emplace_back(as_struct());
  }
  else if (split_ == dim) {
    int_t i = go_right ? 1 : 0;
    child(i)->as_struct()->find_siblings_descend(dim, go_right, siblings);
  }
  else {
    for (auto & c : children_)
      c->as_struct()->find_siblings_descend(dim, go_right, siblings);
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Top level neighbor find
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::find_neighbors()
{

  descend([&](auto & node) {

    const auto & permutations = sector_permutations[num_dims_-1];
    auto struct_node = node.as_struct();

    std::vector<neighbor_pair_t> list;

    for (const auto & sector : permutations) {

      list.clear();
      struct_node->find_neighbors(sector, list);
    
      for (const auto & neigh_pair : list) {
        auto neigh = const_cast<mesh_struct_node_t*>(neigh_pair.first);
        auto root_neigh = neigh_pair.second;
        // different roots
        if (root_neigh) {
        std::cout << root_neigh << std::endl;
          struct_node->add_external_neighbor(
              neigh,
              sector,
              root_neigh->tm_me2donor_cells(),
              root_neigh->tm_me2donor_verts());
        }
        // same roots
        else {
          struct_node->add_internal_neighbor(neigh, sector);
        }
      } // neighbor par

    } // permutations
    
  });
}

////////////////////////////////////////////////////////////////////////////////
/// Top level neighbor find, for a specific sector
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::find_neighbors(
    const sector_t & sector,
    std::vector<mesh_struct_node_t::neighbor_pair_t> & list) const
{
  std::array<bool,3> mask;
  sector.copy(&mask[0]);
  std::cout  << "FIND NEIGHBORS FOR " << sector << " of block " << id_ << std::endl;
  find_neighbors(this, sector, mask, list);
  std::cout  << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
/// Main neighbor finding kernel
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::find_neighbors(
    const mesh_struct_node_t * start,
    const sector_t & sector,
    std::array<bool,3> & mask,
    std::vector<mesh_struct_node_t::neighbor_pair_t> & list) const
{
  auto has_left = std::accumulate(&mask[0], &mask[3], 0);

  // go up until the split matches 
  if (has_left) {

    //--------------------------------------------------------------------------
    // Have a parent
    if (parent_) {

      auto parent = parent_->as_struct();

      // if the split matches, mark it
      auto split_dim = parent->split_;
      auto split_mask = sector[split_dim];
      if (split_mask && mask[split_dim]) {
          
        // if there is a neighbor to cross, cross it
        auto neighbor = (split_mask<0) ? parent->left(this) : parent->right(this);
        if (neighbor) {

          mask[split_dim] = false;
          has_left = std::accumulate(&mask[0], &mask[3], 0);

          if (!has_left) {
            neighbor->find_neighbors_descend(&start->level_[0], start->range_, sector, list);
            return;
          } // none left
          
        } // neighbor

      } // split

      // keep going up if we have more neighbors to find
      parent->find_neighbors(start, sector, mask, list);
    
    } // parent
    //--------------------------------------------------------------------------
    // At a root
    else {

      for (const auto & sector_neighbors : root_neighbors_) {  
        const auto & search_sector = sector_neighbors.first;
        const auto & neighbors = sector_neighbors.second;
        for (const auto & neigh : neighbors) {

          if (!sector.derives_from(search_sector)) continue;

          auto me2donor = neigh.ctm_me2donor_cells();
          me2donor.scale1(&start->level_[0]);

          auto donor_range = start->range_;
          me2donor.transform(donor_range);
          std::cout << "descend " << start->range_ << " vs " << donor_range << std::endl;

          auto donor_sector = sector;
          me2donor.rotate(&donor_sector[0]); 

          auto donor_level = start->level_;
          me2donor.rotate_abs(&donor_level[0]);

          auto end = list.size();
          neigh.node()->as_struct()->find_neighbors_descend(
              &donor_level[0],
              donor_range,
              donor_sector,
              list);

          auto new_end = list.size();
          for (size_t i=end; i<new_end; ++i)
            list[i].second = &neigh;
          std::cout << &neigh << " " << neigh.node() << std::endl;

        } // neighbors
      } // sector neighbors

    } // at a root

  }

}

////////////////////////////////////////////////////////////////////////////////
/// Descend to find neighbors
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::find_neighbors_descend(
    const int_t * level,
    const range_ijk_t & range,
    const sector_t & sector,
    std::vector<mesh_struct_node_t::neighbor_pair_t> & list) const
{
  if (children_.empty()) {
    if (do_ranges_intersect(
          num_dims_,
          sector,
          &level[0],
          range,
          &level_[0],
          range_))
      list.emplace_back(this, nullptr);
  }
  else {
    //--- search dim is same as split, find the closest block to it,
    //--- and stay to one side of the cut plane
    auto split_mask = sector[split_];
    if (split_mask) {

      if (split_mask<0) {
        for (const auto & c : children_) {
          auto chld = c->as_struct();
          if (do_ranges_touch_left(
              split_,
              level,
              &chld->level_[0],
              range,
              chld->range_))
          {
            chld->find_neighbors_descend(level, range, sector, list);
            return;
          }
        }
      } // -1
      else {
        for (const auto & c : children_) {
          auto chld = c->as_struct();
          if (do_ranges_touch_right(
              split_,
              level,
              &chld->level_[0],
              range,
              chld->range_))
          {
            chld->find_neighbors_descend(level, range, sector, list);
            return;
          }
        }
      } // 1

    }
    //--- split along another dim
    else {
      for (const auto & c : children_)
        c->as_struct()->find_neighbors_descend(level, range, sector, list);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Tree refinement
////////////////////////////////////////////////////////////////////////////////
bool mesh_struct_node_t::refine_tree(const amr_flags_t & bflags)
{
  bool did_refine = false;
  for (int_t dim=0; dim<num_dims_; ++dim) {
    if (bflags[dim] == amr_status_t::refine) {
      did_refine = true;
      insert_split(dim);
    }
  }
  return did_refine;
}

////////////////////////////////////////////////////////////////////////////////
/// Tree coarsening
////////////////////////////////////////////////////////////////////////////////
bool mesh_struct_node_t::coarsen_tree(const amr_flags_t & flags)
{
  if (!flags.has_coarsen()) return false;
  
  if (!parent_) {
    std::cout << "CANT COARSEN ROOT" << std::endl;
    return false;
  }

  auto parent = parent_->as_struct();
  auto parent_split = parent->split_;
  
  if (parent_split == -1)
    THROW_ERROR("Parent split -1, something went wrong.");

  auto tmp_flags = flags;
    
  while (tmp_flags.has_coarsen()) {

    // coarsen along the parent split first
    if (flags[parent_split] == amr_status_t::coarsen) {
      tmp_flags[parent_split] = amr_status_t::none;
      assert(
          parent_->child(0)->is_leaf() && parent_->child(1)->is_leaf() &&
          "CANT COARSEN BLOCKS AT DIFFERENT LEVELS");
    }
    // coarsen along other directions
    else {
      auto it = std::find(flags.begin(), flags.end(), amr_status_t::coarsen);
      auto dim = std::distance(flags.begin(), it);
      THROW_ERROR("Not imiplemented yet");
      tmp_flags[dim] = amr_status_t::none;
    }

  } // while

  return true;

}


////////////////////////////////////////////////////////////////////////////////
/// Refine a structured node
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::refine_block(const amr_flags_t & flags)
{
  if (children_.empty()) return;

  // collect ranges
  std::vector<range_ijk_t> ranges;
  std::vector<mesh_struct_node_t*> children;
  
  auto nest = std::pow(2, num_dims_);
  ranges.reserve( nest );
  children.reserve( nest );

  long_t begin[3] = {
    consts::int_max,
    num_dims_>1 ? consts::int_max : 0,
    num_dims_>2 ? consts::int_max : 0};
  
  descend( [&](auto & n){
      auto struct_node = dynamic_cast<mesh_struct_node_t *>(&n);
      const auto & range = struct_node->range_;
      std::cout << "range " << range << std::endl;
      ranges.emplace_back(range); 
      std::cout << "ranges " << ranges.back() << std::endl;
      children.emplace_back(struct_node);
      for (int_t d=0; d<num_dims_; ++d)
        begin[d] = std::min(range.begin[d], begin[d]);
      std::cout  << "dims " << struct_node->dims_[0] << ", " << struct_node->dims_[1] << std::endl;
    } );

  // make ranges relative
  for (auto & r : ranges) {
    r.make_relative( &begin[0] );
    std::cout << "child rage " << r << std::endl;
    for (int_t dim=0; dim<num_dims_; ++dim) {
      auto fact = flags.has_refine(dim) ? 2 : 1;
      r.begin[dim] /= fact;
      r.end[dim] /= fact;
    };
    std::cout << "after rage " << r << std::endl;
  }
  
  // refine
  auto nchilds = ranges.size();
  auto new_blocks = struct_block()->refine(flags, ranges.data(), nchilds);
  auto nnew = new_blocks.size();

  if (nnew != nchilds)
    THROW_ERROR(
        "Tried to refine a block into " << nchilds <<
        ", but returned only " << nnew );

  // move new children to their respective nodes
  for (size_t c=0; c<nchilds; ++c) {
    auto child = children[c];
    child->set_block( std::move(new_blocks[c]) );
  }

}

////////////////////////////////////////////////////////////////////////////////
/// Coarsen a structured node
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_node_t::coarsen_block()
{
  if (children_.empty()) return;

  std::cout << "COARSENING " << std::endl;
  const mesh_struct_block_t* children[] = {
    children_[0]->block()->as_struct(),
    children_[1]->block()->as_struct()
  };
  std::cout << children_[0]->is_leaf() << " "
    << children_[1]->is_leaf() << std::endl;
  
  auto new_block = mesh_struct_block_t::coarsen(split_, children);

  set_block( std::move(new_block) );

}

} // namespace
