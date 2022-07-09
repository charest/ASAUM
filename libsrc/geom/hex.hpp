#ifndef HEX_HPP
#define HEX_HPP

#include "tet.hpp"

namespace prl {
      
namespace detail {

////////////////////////////////////////////////////////////////////////////////
/// Helper to compute the gometry info for a hex[0] cell
////////////////////////////////////////////////////////////////////////////////
inline void _hex_geometry(
    const real_t * ax,
    const real_t * bx,
    const real_t * cx,
    const real_t * dx,
    real_t * xc,
    real_t & v)
{
  // face midpoint
  real_t xm[] = {
    0.25*(ax[0]+bx[0]+cx[0]+dx[0]),
    0.25*(ax[1]+bx[1]+cx[1]+dx[1]),
    0.25*(ax[2]+bx[2]+cx[2]+dx[2])
  };
 
  _tet_geometry(ax, bx, xm, xc, v);
  _tet_geometry(bx, cx, xm, xc, v);
  _tet_geometry(cx, dx, xm, xc, v);
  _tet_geometry(dx, ax, xm, xc, v);
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a hex[0] cell
////////////////////////////////////////////////////////////////////////////////
inline void hex_geometry(
    const real_t * ax,
    const real_t * bx,
    const real_t * cx,
    const real_t * dx,
    const real_t * ex,
    const real_t * fx,
    const real_t * gx,
    const real_t * hx,
    real_t * xc,
    real_t & v)
{
  xc[0] = 0;
  xc[1] = 0;
  xc[2] = 0;
  v = 0;

  // x=i
  detail::_hex_geometry(
    ax, bx, cx,  dx,
    xc, v);
  // x=i+1
  detail::_hex_geometry(
    hx, gx, fx, ex,
    xc, v);
  
  // y=i
  detail::_hex_geometry(
    ax, ex, fx, bx,
    xc, v);
  // y=i+1
  detail::_hex_geometry(
    dx, cx, gx, hx,
    xc, v);
  
  // z=i
  detail::_hex_geometry(
    ax, dx, hx, ex,
    xc, v);
  // z=i+1
  detail::_hex_geometry(
    bx, fx, gx, cx,
    xc, v);
  
  real_t vabs = std::abs(v);
  
  if ( vabs > consts::epsilon ) {
    real_t fact = 1 / (8*v);
    xc[0] *= fact;
    xc[1] *= fact;
    xc[2] *= fact;
  }

  v = vabs / 3;

}

} // prl

#endif
