#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include "config.hpp"
#include <cmath>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a triangular face
////////////////////////////////////////////////////////////////////////////////
inline void triangle_normal(
    const real_t * ax,
    const real_t * bx,
    const real_t * cx,
    real_t * nx)
{
  real_t u[] = { 
    bx[0] - ax[0],
    bx[1] - ax[1],
    bx[2] - ax[2] };
  real_t v[] = {
    cx[0] - ax[0],
    cx[1] - ax[1],
    cx[2] - ax[2] };

  nx[0] = 0.5 * (u[1]*v[2] - u[2]*v[1]);
  nx[1] = 0.5 * (u[2]*v[0] - u[0]*v[2]);
  nx[2] = 0.5 * (u[0]*v[1] - u[1]*v[0]);
}


////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a triangular face
////////////////////////////////////////////////////////////////////////////////
inline void triangle_geometry(
    const real_t * ax,
    const real_t * bx,
    const real_t * cx,
    real_t * fx,
    real_t * nx,
    real_t & a)
{
  constexpr auto third = 1. / 3.;
  fx[0] = third*(ax[0]+bx[0]+cx[0]);
  fx[1] = third*(ax[1]+bx[1]+cx[1]);
  fx[2] = third*(ax[2]+bx[2]+cx[2]);;

  real_t u[] = { 
    bx[0] - ax[0],
    bx[1] - ax[1],
    bx[2] - ax[2] };
  real_t v[] = {
    cx[0] - ax[0],
    cx[1] - ax[1],
    cx[2] - ax[2] };

  nx[0] = 0.5 * (u[1]*v[2] - u[2]*v[1]);
  nx[1] = 0.5 * (u[2]*v[0] - u[0]*v[2]);
  nx[2] = 0.5 * (u[0]*v[1] - u[1]*v[0]);

  a = std::sqrt( nx[0]*nx[0] + nx[1]*nx[1] + nx[2]*nx[2] );
}


} // namesapce

#endif
