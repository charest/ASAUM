#ifndef STRUCTURED_CPU_HPP
#define STRUCTURED_CPU_HPP

#include "config.hpp"

namespace prl{

struct logical_block_t;

/// Update the geometry
void struct_face_geometry(
    const logical_block_t & block,
    const real_t * vx,
    real_t * fx,
    real_t * fn,
    real_t * fa);

/// Update the geometry
void struct_cell_geometry(
    const logical_block_t & block,
    const real_t * vx,
    real_t * xc,
    real_t * cv);
        
// find a cell in a structured block
bool find_cell_struct(
    const logical_block_t & block,
    const real_t * vertices,
    const real_t * x,
    int_t & cell_id);

} // namespace

#endif
