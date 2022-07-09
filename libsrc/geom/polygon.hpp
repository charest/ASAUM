#ifndef POLYGON_HPP
#define POLYGON_HPP

#include "config.hpp"
#include "tri.hpp"

#include <cmath>
#include <iostream>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a polygonal cell
////////////////////////////////////////////////////////////////////////////////
inline void polygon_geometry(
    const real_t * vx,
    const int_t * vs,
    const int_t nv,
    real_t * xc,
    real_t & v)
{
  xc[0] = 0;
  xc[1] = 0;

  v = 0;

  for (int_t i=0; i<nv; ++i) {
    int_t a = vs[i];
    int_t b = vs[(i+1<nv) ? i+1 : 0];
    const real_t * ax = vx + 2*a;
    const real_t * bx = vx + 2*b;

    real_t tmp = ax[0]*bx[1] - bx[0]*ax[1];
    xc[0] += ( ax[0] + bx[0] ) * tmp;
    xc[1] += ( ax[1] + bx[1] ) * tmp;
    v += tmp;
  }  

  real_t vabs = std::abs(v);

  if ( vabs > consts::epsilon ) {
    real_t fact = 1 / (3*v);
    xc[0] *= fact;
    xc[1] *= fact;
  }

  v = 0.5 * vabs;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a polygonal face
////////////////////////////////////////////////////////////////////////////////
inline void polygon_geometry(
    const real_t * vx,
    const int_t * vs,
    const int_t nv,
    real_t * fx,
    real_t * nx,
    real_t & a)
{
  int_t v0;
  real_t fm[] = {0., 0., 0.};
  
  for (int_t i=0; i<nv; ++i) {
    v0 = vs[i];
    fm[0] += vx[ v0*3 ];
    fm[1] += vx[ v0*3 + 1 ];
    fm[2] += vx[ v0*3 + 2 ];
  }
    
  real_t fact = static_cast<real_t>(1) / nv;
  fm[0] *= fact;
  fm[1] *= fact;
  fm[2] *= fact;

  fx[0] = 0;
  fx[1] = 0;
  fx[2] = 0;
  nx[0] = 0;
  nx[1] = 0;
  nx[2] = 0;
  a = 0;

  int_t v1;
  real_t n[3], f[3], atmp;

  for (int_t i=0; i<nv; ++i) {
    v0 = vs[i];
    v1 = (i+1==nv) ? vs[0] : vs[i+1];
    triangle_geometry(
        vx + 3*v0,
        vx + 3*v1,
        fm, f, n, atmp);
    nx[0] += n[0];
    nx[1] += n[1];
    nx[2] += n[2];
    fx[0] += atmp*f[0];
    fx[1] += atmp*f[1];
    fx[2] += atmp*f[2];
    a += atmp;
  }

  if ( a > consts::epsilon ) {
    real_t fact = 1/a;
    fx[0] *= fact;
    fx[1] *= fact;
    fx[2] *= fact;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a polygonal cell
////////////////////////////////////////////////////////////////////////////////
inline bool is_inside_polygon(
    const real_t * vx,
    const int_t * vs,
    const int_t nv,
    const real_t * x)
{

  for (int_t i=0; i<nv; ++i) {
    int_t a = vs[i];
    int_t b = vs[(i+1<nv) ? i+1 : 0];
    const real_t * ax = &vx[2*a];
    const real_t * bx = &vx[2*b];

    real_t dx = bx[0] - ax[0];
    real_t dy = bx[1] - ax[1];

    auto f = dx * (x[1] - ax[1]) - dy * (x[0] - ax[0]);
    if (f<0) return false;
  }  

  return true;

}


} // namesapce

#endif
