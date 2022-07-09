#ifndef UNSTRUCTURED_CPU_HPP
#define UNSTRUCTURED_CPU_HPP

#include "config.hpp"

namespace prl{

/// Update the face geometry
void unstruct_face_geometry_1d(
    int_t nv,
    int_t nc,
    const int_t * cell2vert_offsets,
    const int_t * cell2vert_indices,
    const int_t * vert2cells_offsets,
    const int_t * vert2cells_indices,
    const real_t * vx,
    real_t * fx,
    real_t * fn,
    real_t * fa);
void unstruct_face_geometry_2d(
    int_t nface,
    const int_t * face2vert_offsets,
    const int_t * face2vert_indices,
    const real_t * vx,
    real_t * fx,
    real_t * fn,
    real_t * fa);
void unstruct_face_geometry_3d(
    int_t nface,
    const int_t * face2vert_offsets,
    const int_t * face2vert_indices,
    const real_t * vx,
    real_t * fx,
    real_t * fn,
    real_t * fa);

/// Update the cell geometry
void unstruct_cell_geometry_1d(
    int_t nc,
    const int_t * cell2vert_offsets,
    const int_t * cell2vert_indices,
    const real_t * vx,
    real_t * xc,
    real_t * cv);
void unstruct_cell_geometry_2d(
    int_t nc,
    const int_t * cell2vert_offsets,
    const int_t * cell2vert_indices,
    const real_t * vx,
    real_t * xc,
    real_t * cv);
void unstruct_cell_geometry_3d(
    int_t nc,
    const int_t * cell2faces_offsets,
    const int_t * cell2faces_indices,
    const int_t * face2cells_offsets,
    const int_t * face2cells_indices,
    const int_t * face_owner,
    const real_t * vx,
    real_t * xc,
    real_t * cv);

// find a cell in a unstructured block
bool find_cell_unstruct_1d(
    int_t nc,
    const int_t * cell2verts_offsets,
    const int_t * cell2verts_indices,
    const real_t * vertices,
    const real_t * x,
    int_t & cell_id);
bool find_cell_unstruct_2d(
    int_t nc,
    const int_t * cell2verts_offsets,
    const int_t * cell2verts_indices,
    const real_t * vertices,
    const real_t * x,
    int_t & cell_id);
bool find_cell_unstruct_3d(
    int_t nc,
    const int_t * cell2faces_offsets,
    const int_t * cell2faces_indices,
    const int_t * face2verts_offsets,
    const int_t * face2verts_indices,
    const real_t * vertices,
    const real_t * x,
    int_t & cell_id);

} // namespace

#endif
