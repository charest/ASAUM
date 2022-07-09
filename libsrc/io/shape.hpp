#ifndef SHAPE_HPP
#define SHAPE_HPP

namespace prl {

//! An enumeration to keep track of element types
enum class shape_t {
  line,
  tri,
  quad,
  polygon,
  tet,
  hex,
  polyhedron,
  unknown,
  empty };

} // namespace

#endif
