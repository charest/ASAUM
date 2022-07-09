#ifndef BLOCK_IJK_NEIGHBOR_HPP
#define BLOCK_IJK_NEIGHBOR_HPP

#include "config.hpp"
#include "ctm.hpp"
#include "range.hpp"
#include "sector.hpp"
#include "transform.hpp"

namespace prl {

class block_ijk_t;

////////////////////////////////////////////////////////////////////////////////
/// Neighbor info class 
////////////////////////////////////////////////////////////////////////////////
class block_ijk_neighbor_t {
  const block_ijk_t & neigh_;
  ctm_t me2donor_cells_, me2donor_verts_;
  range_ijk_t range_;

public:

  block_ijk_neighbor_t(
      const block_ijk_t & neigh,
      const transform_t & trans_cells,
      const transform_t & trans_verts,
      const range_ijk_t & my_range);
  
  block_ijk_neighbor_t(
      const block_ijk_t & neigh,
      const ctm_t & trans_cells,
      const ctm_t & trans_verts,
      const range_ijk_t & my_range);

  const block_ijk_t & block() const { return neigh_; }

  const auto & range() const { return range_; }

  transform_t tm_me2donor_cells() const { return me2donor_cells_; }
  transform_t tm_me2donor_verts() const { return me2donor_verts_; }
  
  const ctm_t & ctm_me2donor_cells() const { return me2donor_cells_; }
  const ctm_t & ctm_me2donor_verts() const { return me2donor_verts_; }

  void me2donor_cells(int_t * is) const
  { me2donor_cells_.transform(is); }
  
  void me2donor_verts(int_t * is) const
  { me2donor_verts_.transform(is); }

  bool vert_is_inside(const int_t * pos) const;
};

////////////////////////////////////////////////////////////////////////////////
/// Shared info class 
////////////////////////////////////////////////////////////////////////////////
class block_ijk_shared_t {
  const block_ijk_t & neigh_;
  sector_t sector_;
  range_ijk_t range_;

public:

  block_ijk_shared_t(
      const block_ijk_t & neigh,
      const sector_t & sector,
      const range_ijk_t & my_range) : 
    neigh_(neigh), sector_(sector), range_(my_range)
  {}
  const block_ijk_t & block() const { return neigh_; }

  const auto & sector() const { return sector_; }
  const auto & range() const { return range_; }
};


} // prl

#endif
