#ifndef BLOCK_IJK_VERTEX_T
#define BLOCK_IJK_VERTEX_T

#include "config.hpp"
#include <vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Vertex class for computing block connectivity
////////////////////////////////////////////////////////////////////////////////
struct block_ijk_vertex_t {
  int_t id_ = -1;
  std::vector<std::vector<real_t>> dirs_;
  
  block_ijk_vertex_t(int_t id, int_t ndims)
    : id_(id), dirs_(ndims, std::vector<real_t>(ndims))
  {}

  real_t * direction(int_t d) { return dirs_[d].data(); }
};


} // namespace


#endif
