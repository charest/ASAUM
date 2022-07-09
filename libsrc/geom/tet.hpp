#ifndef TET_HPP
#define TET_HPP

#include "tri.hpp"
#include "math/core.hpp"

namespace prl {

namespace detail {

////////////////////////////////////////////////////////////////////////////////
/// Helper to compute the gometry info for a hex[0] cell
////////////////////////////////////////////////////////////////////////////////
inline void _tet_geometry(
    const real_t * ax,
    const real_t * bx,
    const real_t * cx,
    real_t * xc,
    real_t & v)
{
  real_t n[3];
  triangle_normal(ax, bx, cx, n);
  xc[0] += n[0] * (sqr(ax[0] + bx[0]) + sqr(bx[0] + cx[0]) + sqr(ax[0] + cx[0]));
  xc[1] += n[1] * (sqr(ax[1] + bx[1]) + sqr(bx[1] + cx[1]) + sqr(ax[1] + cx[1]));
  xc[2] += n[2] * (sqr(ax[2] + bx[2]) + sqr(bx[2] + cx[2]) + sqr(ax[2] + cx[2]));
  v += n[0]*cx[0] + n[1]*cx[1] + n[2]*cx[2];
}

} // namespace

} // namespace

#endif
