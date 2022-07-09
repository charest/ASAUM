#ifndef RANGE_HPP
#define RANGE_HPP

#include "config.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <ostream>

namespace prl {

struct range_ijk_t;
class sector_t;

bool do_ranges_intersect(
    const int_t nd,
    const sector_t & sector,
    const range_ijk_t & my_range,
    const range_ijk_t & neigh_range);

bool intersect_ranges(
    const int_t nd,
    const sector_t & sector,
    range_ijk_t & my_range,
    range_ijk_t & neigh_range);

////////////////////////////////////////////////////////////////////////////////
/// ijk range class
////////////////////////////////////////////////////////////////////////////////
struct range_ijk_t
{
  std::array<long_t,3> begin = {0, 0, 0};
  std::array<long_t,3> end = {1, 1, 1};
  
  range_ijk_t() = default;

  range_ijk_t(int_t nd, const long_t * b, const long_t * e)
  {
    std::copy_n(b, nd, &begin[0]);
    std::copy_n(e, nd, &end[0]);
  }
  
  range_ijk_t(int_t nd, const long_t * dims)
  { std::copy_n(dims, nd, &end[0]); }
  
  auto size(int_t i) const { return end[i] - begin[i]; }

  template<typename T>
  void mult(int_t i, T x)
  {
    begin[i] *= x;
    end[i] *= x;
  }

  template<typename T>
  void div(int_t i, T x)
  {
    begin[i] /= x;
    end[i] /= x;
  }

  friend std::ostream &operator<<( std::ostream &out, const range_ijk_t & r )
  {
    out << "[(" << r.begin[0]
        << ", " << r.begin[1]
        << ", " << r.begin[2] << "), "
        << "("  << r.end[0]
        << ", " << r.end[1]
        << ", " << r.end[2] << ")]"; 
    return out;
  }

  /// Is an index inside a range
  template<typename T>
  bool is_inside(const T * vals, size_t n, bool include_end = false) const
  {
    if (include_end) {
      for (size_t i=0; i<n; ++i)
        if (vals[i]<begin[i] || vals[i]>end[i])
          return false;
    }
    else {
      for (size_t i=0; i<n; ++i)
        if (vals[i]<begin[i] || vals[i]>=end[i])
          return false;
    }
    return true;
  }

  void make_relative(long_t * pos);
};

} // namespace

#endif
