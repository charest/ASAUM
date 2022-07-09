#include "config.hpp"
#include "geom/polygon.hpp"
#include "geom/polyhedron.hpp"
#include "math/core.hpp"

#include <cassert>
#include <cmath>

#include <iostream>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Update the face geometry (1d)
////////////////////////////////////////////////////////////////////////////////
void unstruct_face_geometry_1d(
    int_t nv,
    int_t nc,
    const int_t * cell2verts_offsets,
    const int_t * cell2verts_indices,
    const int_t * vert2cells_offsets,
    const int_t * vert2cells_indices,
    const real_t * vx,
    real_t * fx,
    real_t * fn,
    real_t * fa)
{
  
  for (int_t v=0; v<nv; ++v) {
    auto cell_start = vert2cells_offsets[v];
    auto cell_end = vert2cells_offsets[v+1];
    auto n = cell_end - cell_start;
    // physcial/processor boundary
    if (n==1 || (n==2 && vert2cells_indices[cell_start+1] >= nc)) {
      auto c = vert2cells_indices[cell_start];
      auto vert_start = cell2verts_offsets[c];
      auto v0 = cell2verts_indices[vert_start];
      auto v1 = cell2verts_indices[vert_start+1];
      if (v0==v)
        fn[v] = sgn(vx[v0] - vx[v1]);
      else
        fn[v] = sgn(vx[v1] - vx[v0]);
    }
    // internal
    else {
      auto c0 = vert2cells_indices[cell_start];
      auto c1 = vert2cells_indices[cell_start+1];
      auto vert_start = cell2verts_offsets[c0];
      auto v0 = cell2verts_indices[vert_start];
      auto v1 = cell2verts_indices[vert_start+1];
      auto vl = (v0==v) ? v1 : v0;
      vert_start = cell2verts_offsets[c1];
      v0 = cell2verts_indices[vert_start];
      v1 = cell2verts_indices[vert_start+1];
      auto vr = (v0==v) ? v1 : v0;
      fn[v] = sgn(vx[vr] - vx[vl]);
    }
    fx[v] = vx[v];
    fa[v] = 1;
  }

}

////////////////////////////////////////////////////////////////////////////////
/// Update the unstructured face geometry (2d)
////////////////////////////////////////////////////////////////////////////////
void unstruct_face_geometry_2d(
    int_t nface,
    const int_t * face2vert_offsets,
    const int_t * face2vert_indices,
    const real_t * vx,
    real_t * fx,
    real_t * fn,
    real_t * fa)
{
  for (int_t f=0; f<nface; ++f) {
    auto vert_start = face2vert_offsets[f];
    assert(face2vert_offsets[f+1]-vert_start == 2 &&
        "Face should have 2 verts");
    auto a = face2vert_indices[vert_start];
    auto b = face2vert_indices[vert_start+1];
    real_t ax = vx[a*2], ay = vx[a*2+1];
    real_t bx = vx[b*2], by = vx[b*2+1];
    fx[f*2 + 0] = 0.5 * (ax + bx);
    fx[f*2 + 1] = 0.5 * (ay + by);
    auto nx = by - ay;
    auto ny = ax - bx;
    fn[f*2 + 0] = nx;
    fn[f*2 + 1] = ny;
    fa[f] = std::sqrt( nx*nx + ny*ny );
    if (fa[f] > consts::epsilon) { 
      fn[f*2 + 0] /= fa[f];
      fn[f*2 + 1] /= fa[f];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Update the unstructured face geometry (3d)
////////////////////////////////////////////////////////////////////////////////
void unstruct_face_geometry_3d(
    int_t nface,
    const int_t * face2vert_offsets,
    const int_t * face2vert_indices,
    const real_t * vx,
    real_t * fx,
    real_t * fn,
    real_t * fa)
{
  for (int_t f=0; f<nface; ++f) {
    auto vert_start = face2vert_offsets[f];
    auto vert_end = face2vert_offsets[f+1];
    auto nvert = vert_end - vert_start;
    assert(nvert > 2 && "Face should have more than 2 verts");
    polygon_geometry(
        vx,
        face2vert_indices + vert_start,
        nvert,
        fx + f*3,
        fn + f*3,
        fa[f]);
    if (fa[f] > consts::epsilon) { 
      fn[f*3 + 0] /= fa[f];
      fn[f*3 + 1] /= fa[f];
      fn[f*3 + 2] /= fa[f];
    }
  } // face
}

////////////////////////////////////////////////////////////////////////////////
/// Update the cell geometry (1d)
////////////////////////////////////////////////////////////////////////////////
void unstruct_cell_geometry_1d(
    int_t nc,
    const int_t * cell2vert_offsets,
    const int_t * cell2vert_indices,
    const real_t * vx,
    real_t * xc,
    real_t * cv)
{

  for (int_t c=0; c<nc; ++c) {
    auto start = cell2vert_offsets[c];
    assert(cell2vert_offsets[c+1]-start==2 && "Cell should have 2 verts");
    auto v0 = cell2vert_indices[start];
    auto v1 = cell2vert_indices[start+1];
    xc[c] = 0.5 * (vx[v1] + vx[v0]);
    cv[c] = std::abs(vx[v1] - vx[v0]);
  }
  
}

////////////////////////////////////////////////////////////////////////////////
/// Update the geometry (2d)
////////////////////////////////////////////////////////////////////////////////
void unstruct_cell_geometry_2d(
    int_t nc,
    const int_t * cell2vert_offsets,
    const int_t * cell2vert_indices,
    const real_t * vx,
    real_t * xc,
    real_t * cv)
{

  for (int_t c=0; c<nc; ++c) {
    auto start = cell2vert_offsets[c];
    auto end = cell2vert_offsets[c+1];
    auto nv = end - start;
    assert(nv>2 && "Cell should more than 2 verts");
    polygon_geometry(
        vx,
        &cell2vert_indices[start],
        nv,
        xc + c*2,
        cv[c]);
  }
  
}

////////////////////////////////////////////////////////////////////////////////
/// Update the geometry (3d)
////////////////////////////////////////////////////////////////////////////////
void unstruct_cell_geometry_3d(
    int_t nc,
    const int_t * cell2face_offsets,
    const int_t * cell2face_indices,
    const int_t * face2vert_offsets,
    const int_t * face2vert_indices,
    const int_t * face_owner,
    const real_t * vx,
    real_t * xc,
    real_t * cv)
{

  for (int_t c=0; c<nc; ++c) {
    auto start = cell2face_offsets[c];
    auto end = cell2face_offsets[c+1];
    auto nf = end - start;
    assert(nf>2 && "Cell should more than 2 faces");
    polyhedron_geometry(
        vx,
        face2vert_offsets,
        face2vert_indices,
        c,
        face_owner,
        &cell2face_indices[start],
        nf,
        xc + c*3,
        cv[c]);
  }
  
}

////////////////////////////////////////////////////////////////////////////////
// find a cell in a unstructured block
////////////////////////////////////////////////////////////////////////////////
bool find_cell_unstruct_1d(
    int_t nc,
    const int_t * cell2verts_offsets,
    const int_t * cell2verts_indices,
    const real_t * vertices,
    const real_t * x,
    int_t & cell_id)
{
  for (cell_id=0; cell_id<nc; ++cell_id) {
    auto start = cell2verts_offsets[cell_id];
    auto sz = cell2verts_offsets[cell_id-1] - start;
    if (sz != 2) continue;
    auto a = cell2verts_indices[start];
    auto b = cell2verts_indices[start+1];
    if (!(x[0] < vertices[a] && x[0] > vertices[b]) ||
        !(x[0] < vertices[b] && x[0] > vertices[a]))
      return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// find a cell in a unstructured block
////////////////////////////////////////////////////////////////////////////////
bool find_cell_unstruct_2d(
    int_t nc,
    const int_t * cell2verts_offsets,
    const int_t * cell2verts_indices,
    const real_t * vertices,
    const real_t * x,
    int_t & cell_id)
{
  for (cell_id=0; cell_id<nc; ++cell_id) {
    auto start = cell2verts_offsets[cell_id];
    auto sz = cell2verts_offsets[cell_id-1] - start;
    if (sz < 3) continue;
    auto is_inside = is_inside_polygon(vertices, &cell2verts_indices[start], sz, x);
    if (is_inside) return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// find a cell in a unstructured block
////////////////////////////////////////////////////////////////////////////////
bool find_cell_unstruct_3d(
    int_t nc,
    const int_t * cell2faces_offsets,
    const int_t * cell2faces_indices,
    const int_t * face2verts_offsets,
    const int_t * face2verts_indices,
    const real_t * vertices,
    const real_t * x,
    int_t & cell_id)
{
  for (cell_id=0; cell_id<nc; ++cell_id) {
    auto start = cell2faces_offsets[cell_id];
    auto nf = cell2faces_offsets[cell_id-1] - start;
    if (nf < 3) continue;
    auto is_inside = is_inside_polyhedron(
        vertices,
        &cell2faces_indices[start],
        face2verts_offsets,
        face2verts_indices,
        nf,
        x);
    if (is_inside) return true;
  }
  return false;
}

} // namespace
