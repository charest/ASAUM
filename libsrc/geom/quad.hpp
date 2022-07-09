#ifndef QUAD_HPP
#define QUAD_HPP

#include "tri.hpp"

#include <cmath>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a quad face
////////////////////////////////////////////////////////////////////////////////
inline void quad_geometry(
    const real_t * ax,
    const real_t * bx,
    const real_t * cx,
    const real_t * dx,
    real_t * fx,
    real_t * nx,
    real_t & a)
{
  real_t fm[] = {
    0.25*(ax[0] + bx[0] + cx[0] + dx[0]),
    0.25*(ax[1] + bx[1] + cx[1] + dx[1]),
    0.25*(ax[2] + bx[2] + cx[2] + dx[2])};

  fx[0] = 0;
  fx[1] = 0;
  fx[2] = 0;
  nx[0] = 0;
  nx[1] = 0;
  nx[2] = 0;
  a = 0;

  real_t n[3], f[3], atmp;

  triangle_geometry(ax, bx, fm, f, n, atmp);
  nx[0] += n[0];
  nx[1] += n[1];
  nx[2] += n[2];
  fx[0] += atmp*f[0];
  fx[1] += atmp*f[1];
  fx[2] += atmp*f[2];
  a += atmp;
  
  triangle_geometry(bx, cx, fm, f, n, atmp);
  nx[0] += n[0];
  nx[1] += n[1];
  nx[2] += n[2];
  fx[0] += atmp*f[0];
  fx[1] += atmp*f[1];
  fx[2] += atmp*f[2];
  a += atmp;
  
  triangle_geometry(cx, dx, fm, f, n, atmp);
  nx[0] += n[0];
  nx[1] += n[1];
  nx[2] += n[2];
  fx[0] += atmp*f[0];
  fx[1] += atmp*f[1];
  fx[2] += atmp*f[2];
  a += atmp;
  
  triangle_geometry(dx, ax, fm, f, n, atmp);
  nx[0] += n[0];
  nx[1] += n[1];
  nx[2] += n[2];
  fx[0] += atmp*f[0];
  fx[1] += atmp*f[1];
  fx[2] += atmp*f[2];
  a += atmp;

  if ( a > consts::epsilon ) {
    real_t fact = 1/a;
    fx[0] *= fact;
    fx[1] *= fact;
    fx[2] *= fact;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a quad cell
////////////////////////////////////////////////////////////////////////////////
inline void quad_geometry(
    const real_t * ax,
    const real_t * bx,
    const real_t * cx,
    const real_t * dx,
    real_t * xc,
    real_t & v)
{
  xc[0] = 0;
  xc[1] = 0;

  v = 0;

  auto tmp = ax[0]*bx[1] - bx[0]*ax[1];
  xc[0] += ( ax[0] + bx[0] ) * tmp;
  xc[1] += ( ax[1] + bx[1] ) * tmp;
  v += tmp;
  
  tmp = bx[0]*cx[1] - cx[0]*bx[1];
  xc[0] += ( bx[0] + cx[0] ) * tmp;
  xc[1] += ( bx[1] + cx[1] ) * tmp;
  v += tmp;
  
  tmp = cx[0]*dx[1] - dx[0]*cx[1];
  xc[0] += ( cx[0] + dx[0] ) * tmp;
  xc[1] += ( cx[1] + dx[1] ) * tmp;
  v += tmp;
  
  tmp = dx[0]*ax[1] - ax[0]*dx[1];
  xc[0] += ( dx[0] + ax[0] ) * tmp;
  xc[1] += ( dx[1] + ax[1] ) * tmp;
  v += tmp;

  real_t vabs = std::abs(v);

  if ( vabs > consts::epsilon ) {
    real_t fact = 1 / (3*v);
    xc[0] *= fact;
    xc[1] *= fact;
  }

  v = 0.5 * vabs;
}


} // namesapce

#endif
