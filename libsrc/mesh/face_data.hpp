#ifndef INTRABLOCK_FACE_HPP
#define INTRABLOCK_FACE_HPP

#include "config.hpp"

namespace prl {

struct boundary_face_t {
  int_t id;
  int_t dir;
  int_t owner;
  boundary_face_t(int_t i, int_t d, int_t o) : id(i), dir(d), owner(o) {}
};


struct intrablock_face_t {
  int_t id;
  int_t dir;
  int_t owner, ghost;
  intrablock_face_t(int_t i, int_t d, int_t own, int_t gho) 
    : id(i), dir(d), owner(own), ghost(gho) {}
};

} // namespace

#endif
