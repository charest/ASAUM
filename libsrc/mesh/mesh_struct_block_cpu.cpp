#include "config.hpp"
#include "logical/logical_iterators.hpp"

#include "geom/quad.hpp"
#include "geom/hex.hpp"
#include "geom/polygon.hpp"
#include "geom/polyhedron.hpp"


namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Update the geometry (3d)
////////////////////////////////////////////////////////////////////////////////
void struct_face_geometry(
    const logical_block_t & block,
    const real_t * vx,
    real_t * fx,
    real_t * fn,
    real_t * fa)
{

  auto nd = block.nd;
  auto nx = block.ni;
  auto ny = block.nj;
  auto nz = block.nk;

  const auto vertices = block.vertices();
  const logical_face_t faces[] = {
    block.faces(0),
    block.faces(1),
    block.faces(2)
  };

  //------------------------------------
  if (nd == 1) {
  
    auto dir = sgn(vx[1] - vx[0]);
    for (int_t i=0; i<nx+1; ++i) {
      fx[i] = vx[i];
      fn[i] = dir;
      fa[i] = 1;
    }

  }
  
  //------------------------------------
  else if (nd == 2) {

    auto face = &faces[0];

    for (int_t j=0; j<ny; ++j) {
      for (int_t i=0; i<nx+1; ++i) {
        auto fid = face->id(i, j, 0);
        auto a = vertices.id(i, j  , 0);
        auto b = vertices.id(i, j+1, 0);
        real_t ax = vx[a*2], ay = vx[a*2+1];
        real_t bx = vx[b*2], by = vx[b*2+1];
        fx[fid*2 + 0] = 0.5 * (ax + bx);
        fx[fid*2 + 1] = 0.5 * (ay + by);
        auto nx = by - ay;
        auto ny = ax - bx;
        fn[fid*2 + 0] = nx;
        fn[fid*2 + 1] = ny;
        fa[fid] = std::sqrt( nx*nx + ny*ny );
        if (fa[fid] > consts::epsilon) { 
          fn[fid*2 + 0] /= fa[fid];
          fn[fid*2 + 1] /= fa[fid];
        }
      }
    }

    face = &faces[1];
    
    for (int_t j=0; j<ny+1; ++j) {
      for (int_t i=0; i<nx; ++i) {
        auto fid = face->id(i, j, 0);
        auto a = vertices.id(i  , j, 0);
        auto b = vertices.id(i+1, j, 0);
        real_t ax = vx[a*2], ay = vx[a*2+1];
        real_t bx = vx[b*2], by = vx[b*2+1];
        fx[fid*2 + 0] = 0.5 * (ax + bx);
        fx[fid*2 + 1] = 0.5 * (ay + by);
        auto nx = ay - by;
        auto ny = bx - ax;
        fn[fid*2 + 0] = nx;
        fn[fid*2 + 1] = ny;
        fa[fid] = std::sqrt( nx*nx + ny*ny );
        if (fa[fid] > consts::epsilon) { 
          fn[fid*2 + 0] /= fa[fid];
          fn[fid*2 + 1] /= fa[fid];
        }
      }
    }
  }
  //------------------------------------
  else {
    
    auto face = &faces[0];

    for (int_t k=0; k<nz; ++k) {
      for (int_t j=0; j<ny; ++j) {
        for (int_t i=0; i<nx+1; ++i) {
          auto fid = face->id(i, j, k);
          auto a = vertices.id(i, j+0, k+0);
          auto b = vertices.id(i, j+1, k+0);
          auto c = vertices.id(i, j+1, k+1);
          auto d = vertices.id(i, j+0, k+1);
          quad_geometry(
              vx + a*3,
              vx + b*3,
              vx + c*3,
              vx + d*3,
              fx + fid*3,
              fn + fid*3,
              fa[fid]);
          if (fa[fid] > consts::epsilon) { 
            fn[fid*3 + 0] /= fa[fid];
            fn[fid*3 + 1] /= fa[fid];
            fn[fid*3 + 2] /= fa[fid];
          }
        }
      }
    }
    
    face = &faces[1];
    
    for (int_t k=0; k<nz; ++k) {
      for (int_t j=0; j<ny+1; ++j) {
        for (int_t i=0; i<nx; ++i) {
          auto fid = face->id(i, j, k);
          auto a = vertices.id(i+0, j, k+0);
          auto b = vertices.id(i+0, j, k+1);
          auto c = vertices.id(i+1, j, k+1);
          auto d = vertices.id(i+1, j, k+0);
          quad_geometry(
              vx + a*3,
              vx + b*3,
              vx + c*3,
              vx + d*3,
              fx + fid*3,
              fn + fid*3,
              fa[fid]);
          if (fa[fid] > consts::epsilon) { 
            fn[fid*3 + 0] /= fa[fid];
            fn[fid*3 + 1] /= fa[fid];
            fn[fid*3 + 2] /= fa[fid];
          }
        }
      }
    }
    
    face = &faces[2];
    
    for (int_t k=0; k<nz+1; ++k) {
      for (int_t j=0; j<ny; ++j) {
        for (int_t i=0; i<nx; ++i) {
          auto fid = face->id(i, j, k);
          auto a = vertices.id(i+0, j+0, k);
          auto b = vertices.id(i+1, j+0, k);
          auto c = vertices.id(i+1, j+1, k);
          auto d = vertices.id(i+0, j+1, k);
          quad_geometry(
              vx + a*3,
              vx + b*3,
              vx + c*3,
              vx + d*3,
              fx + fid*3,
              fn + fid*3,
              fa[fid]);
          if (fa[fid] > consts::epsilon) { 
            fn[fid*3 + 0] /= fa[fid];
            fn[fid*3 + 1] /= fa[fid];
            fn[fid*3 + 2] /= fa[fid];
          }
        }
      }
    }

  } // num_dims

}

////////////////////////////////////////////////////////////////////////////////
/// Update the geometry (3d)
////////////////////////////////////////////////////////////////////////////////
void struct_cell_geometry(
    const logical_block_t & block,
    const real_t * vx,
    real_t * xc,
    real_t * cv)
{
  auto nd = block.nd;
  auto vertices = block.vertices();
  auto cells = block.cells();
  auto cit = cells.begin();

  int_t vs[8] = {0, 0, 0, 0, 0, 0, 0, 0};

  //------------------------------------
  if (nd == 1) {

    do {
      cit.vertex_ids(vertices, vs);
      auto id = cit.id();
      auto a = vs[0];
      auto b = vs[1];
      xc[id] = 0.5 * (vx[b] + vx[a]);
      cv[id] = std::abs(vx[b] - vx[a]);
    } while(cit.next());

  }
  //------------------------------------
  else if (nd == 2) {
    
    do {
      cit.vertex_ids(vertices, vs);
      auto id = cit.id();
      quad_geometry(
          vx + vs[0]*2,
          vx + vs[3]*2,
          vx + vs[2]*2,
          vx + vs[1]*2,
          xc + id*2,
          cv[id]);
    } while(cit.next());

  }
  //------------------------------------
  else {
    
    do {
      cit.vertex_ids(vertices, vs);
      auto id = cit.id();
      hex_geometry(
          vx + vs[0]*3,
          vx + vs[3]*3,
          vx + vs[2]*3,
          vx + vs[1]*3,
          vx + vs[4]*3,
          vx + vs[7]*3,
          vx + vs[6]*3,
          vx + vs[5]*3,
          xc + id*3,
          cv[id]);
    } while(cit.next());

  } // num_dims
  
}


////////////////////////////////////////////////////////////////////////////////
// find a cell in a structured block
////////////////////////////////////////////////////////////////////////////////
bool find_cell_struct(
    const logical_block_t & block,
    const real_t * coords,
    const real_t * x,
    int_t & cell_id)
{

  auto nd = block.nd;

  auto cells = block.cells();
  auto vertices = block.vertices();

  auto cit = cells.begin();
  
  do {

    cell_id = cit.id();
    
    //--- one dimensional
    if (nd==1) {
      auto a = cit.vertex_id(vertices, 0, 0, 0);
      auto b = cit.vertex_id(vertices, 1, 0, 0);
      if (!(x[0] < coords[a] && x[0] > coords[b]) ||
          !(x[0] < coords[b] && x[0] > coords[a]))
        return true;
    }
    //--- two dimensional
    else if (nd==2) {
      int_t cell_vs[4] = {
        cit.vertex_id(vertices, 0, 0, 0),
        cit.vertex_id(vertices, 1, 0, 0),
        cit.vertex_id(vertices, 1, 1, 0),
        cit.vertex_id(vertices, 0, 1, 0) 
      };
      auto is_inside = is_inside_polygon(coords, cell_vs, 4, x);
      if (is_inside) return true;
    }
    //--- three dimensional
    else {
      auto a = cit.vertex_id(vertices, 0, 0, 0);
      auto b = cit.vertex_id(vertices, 1, 0, 0);
      auto c = cit.vertex_id(vertices, 1, 1, 0);
      auto d = cit.vertex_id(vertices, 0, 1, 0);
      auto e = cit.vertex_id(vertices, 0, 0, 1);
      auto f = cit.vertex_id(vertices, 1, 0, 1);
      auto g = cit.vertex_id(vertices, 1, 1, 1);
      auto h = cit.vertex_id(vertices, 0, 1, 1);
      int_t cell_faces[6] = {0, 1, 2, 3, 4, 5};      
      int_t face2vert_offsets[7] = {0, 4, 8, 12, 16, 20, 24};
      int_t face2vert_indices[] = {
        a, b, c, d,
        a, e, f, b,
        b, f, g, c,
        c, g, h, d,
        d, h, e, a,
        e, h, g, f
      };
      auto is_inside = is_inside_polyhedron(
          coords,
          cell_faces,
          face2vert_offsets,
          face2vert_indices,
          6,
          x);
      if (is_inside) return true;
    }


  } while(cit.next());

  return false;
}

} // namespace
