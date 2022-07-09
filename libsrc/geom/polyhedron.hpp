#ifndef POLYHRDRON_HPP
#define POLYHRDRON_HPP

#include "config.hpp"
#include "hex.hpp"
#include "tri.hpp"
#include "math/math.hpp"

namespace prl {
      
namespace detail {

////////////////////////////////////////////////////////////////////////////////
/// Helper to compute the gometry info for a hex[0] cell
////////////////////////////////////////////////////////////////////////////////
inline void _polyhedron_geometry(
    const real_t * vx,
    const int_t * vs,
    int_t nv,
    real_t * xc,
    real_t & v,
    bool flip)
{
  // face midpoint
  int_t v0;
  real_t xm[] = {0., 0., 0.};
  
  for (int_t i=0; i<nv; ++i) {
    v0 = vs[i];
    xm[0] += vx[ v0*3 ];
    xm[1] += vx[ v0*3 + 1 ];
    xm[2] += vx[ v0*3 + 2 ];
  }
    
  real_t fact = static_cast<real_t>(1) / nv;
  xm[0] *= fact;
  xm[1] *= fact;
  xm[2] *= fact;
  
  int_t v1;

  if (flip) {
    for (int_t i=nv; i --> 0; ) {
      v0 = vs[i];
      v1 = (i==0) ? vs[nv-1] : vs[i-1];
      _tet_geometry(
          vx + 3*v0,
          vx + 3*v1,
          xm, xc, v);
    }
  }
  else {
    for (int_t i=0; i<nv; ++i) {
      v0 = vs[i];
      v1 = (i+1==nv) ? vs[0] : vs[i+1];
      _tet_geometry(
          vx + 3*v0,
          vx + 3*v1,
          xm, xc, v);
    }
  }
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
/// Compute the gometry info for a hex[0] cell
////////////////////////////////////////////////////////////////////////////////
inline void polyhedron_geometry(
    const real_t * vx,
    const int_t * f2v_offsets,
    const int_t * f2v_indices,
    int_t c,
    const int_t * face_owner,
    const int_t * fs,
    int_t nf,
    real_t * xc,
    real_t & v)
{
  xc[0] = 0;
  xc[1] = 0;
  xc[2] = 0;
  v = 0;

  bool flip;
  int_t start, end, nv, fid;
    
  for (int_t f=0; f<nf; ++f) {
    fid = fs[f];
    start = f2v_offsets[fid];
    end = f2v_offsets[fid+1];
    nv = end - start;
    flip = face_owner[fid] != c;
    detail::_polyhedron_geometry(
      vx, f2v_indices + start,
      nv, xc, v, flip);
  }

  real_t vabs = std::abs(v);
  
  if ( vabs > consts::epsilon ) {
    real_t fact = 1 / (8*v);
    xc[0] *= fact;
    xc[1] *= fact;
    xc[2] *= fact;
  }

  v = vabs / 3;

}

////////////////////////////////////////////////////////////////////////////////
/// Is a point inside a polyhedron
////////////////////////////////////////////////////////////////////////////////
inline bool is_inside_polyhedron(
    const real_t * vx,
    const int_t * fs,
    const int_t * f2v_offsets,
    const int_t * f2v_indices,
    int_t nf,
    const real_t * x)
{
  real_t fx[3];
  real_t nx[3];
  real_t dx[3];

  for (int_t f=0; f<nf; ++f) {
  
    auto fid = fs[f];
    auto start = f2v_offsets[fid];
    auto nv = f2v_offsets[fid+1] - start;

    //--- midpoint
    fx[0] = 0.;
    fx[1] = 0.;
    fx[2] = 0.;

    for (auto v=0; v<nv; ++v) {
      auto vid = f2v_indices[start + v];
      auto xv = &vx[vid*3];
      for (int_t d=0; d<3; ++d)
        fx[d] += xv[d];
    }
    for (int_t d=0; d<3; ++d)
      fx[d] /= nv;
      
    dx[0] = x[0] - fx[0];
    dx[1] = x[1] - fx[1];
    dx[2] = x[2] - fx[2];

    //--- side
    for (auto v0=0; v0<nv; ++v0) {
      auto v1 = v0+1<nv ? v0+1 : 0;
      auto vid0 = f2v_indices[start + v0];
      auto vid1 = f2v_indices[start + v1];
      auto xv0 = &vx[vid0*3];
      auto xv1 = &vx[vid1*3];
      triangle_normal(xv0, xv1, fx, nx);
      auto f = dot_nd(nx, dx, 3);
      if (f<0) return false;
    }
  }

  // success if got here
  return true;

}

} // prl

#endif
