#include "block_ijk_neighbor.hpp"
#include "block_ijk.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Neighbor class cconstructor
////////////////////////////////////////////////////////////////////////////////
block_ijk_neighbor_t::block_ijk_neighbor_t(
    const block_ijk_t & neigh,
    const transform_t & trans_cells,
    const transform_t & trans_verts,
    const range_ijk_t & my_range) :
  neigh_(neigh),
  me2donor_cells_(trans_cells.num_dims()),
  me2donor_verts_(trans_verts.num_dims()),
  range_(my_range)
{
  trans_cells.compact(me2donor_cells_);
  trans_verts.compact(me2donor_verts_);
}

block_ijk_neighbor_t::block_ijk_neighbor_t(
    const block_ijk_t & neigh,
    const ctm_t & trans_cells,
    const ctm_t & trans_verts,
    const range_ijk_t & my_range) :
  neigh_(neigh),
  me2donor_cells_(trans_cells),
  me2donor_verts_(trans_verts),
  range_(my_range)
{}

////////////////////////////////////////////////////////////////////////////////
/// Check if neighbor is inside
////////////////////////////////////////////////////////////////////////////////
bool block_ijk_neighbor_t::vert_is_inside(const int_t * pos) const
{ return range_.is_inside(pos, neigh_.num_dims(), true); }

} // namespace
