#include "range.hpp"
#include "sector.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Intersect two bounding boxes (sort of)
////////////////////////////////////////////////////////////////////////////////
bool do_ranges_intersect(
    const int_t nd,
    const sector_t & sector,
    const range_ijk_t & my_range,
    const range_ijk_t & neigh_range)
{
  for (int_t dim=0; dim<nd; ++dim) {
    if (sector[dim] == 0) {
      auto beg = std::max(my_range.begin[dim], neigh_range.begin[dim]);
      auto end = std::min(my_range.end[dim], neigh_range.end[dim]);
      auto len = end - beg;
      if (len<1) return false;
    }
    else {
      if (sector[dim] < 0) {
        if (neigh_range.begin[dim] >= my_range.begin[dim] ||
            neigh_range.end[dim]   <  my_range.begin[dim])
          return false;
      }
      else /*if (sector[dim] > 0)*/ {
        if (neigh_range.begin[dim] >  my_range.end[dim] ||
            neigh_range.end[dim]   <= my_range.end[dim])
          return false;
      }
    } // else
  } // dim

  return true;
}


////////////////////////////////////////////////////////////////////////////////
/// Intersect two bounding boxes (sort of)
////////////////////////////////////////////////////////////////////////////////
bool intersect_ranges(
    const int_t nd,
    const sector_t & sector,
    range_ijk_t & my_range,
    range_ijk_t & neigh_range)
{
  for (int_t dim=0; dim<nd; ++dim) {
    if (sector[dim] == 0) {
      auto beg = std::max(my_range.begin[dim], neigh_range.begin[dim]);
      auto end = std::min(my_range.end[dim], neigh_range.end[dim]);
      neigh_range.begin[dim] = beg;
      my_range   .begin[dim] = beg;
      neigh_range.end  [dim] = end;
      my_range   .end  [dim] = end;
      auto len = end - beg;
      if (len<1) return false;
    }
    else {
      if (sector[dim] < 0) {
        if (neigh_range.begin[dim] >= my_range.begin[dim] ||
            neigh_range.end[dim]   <  my_range.begin[dim])
          return false;
        neigh_range.end[dim] = std::min( my_range.begin[dim], neigh_range.end[dim] );
      }
      else if (sector[dim] > 0) {
        if (neigh_range.begin[dim] >  my_range.end[dim] ||
            neigh_range.end[dim]   <= my_range.end[dim])
          return false;
        neigh_range.begin[dim] = std::max( my_range.end[dim], neigh_range.begin[dim]);
      }
    } // else
  } // dim

  return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Make a range relative to me
////////////////////////////////////////////////////////////////////////////////
void range_ijk_t::make_relative(long_t * pos)
{
  for (int_t dim=0; dim<3; ++dim) {
    begin[dim] -= pos[dim];
    end[dim] -= pos[dim];
  }
}


} // namespace
