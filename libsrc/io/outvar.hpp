#ifndef OUTVAR_HPP
#define OUTVAR_HPP

#include "config.hpp"

#include <string>

namespace prl {

struct outvar_t {
  std::string name;
  int components;
  bool dimensional;
  outvar_t(const std::string & lbl, int comp=1, bool is_dim=false) :
    name(lbl), components(comp), dimensional(is_dim)
  {}
};


} // namespace

#endif
