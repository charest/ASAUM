#ifndef BAMR_FLAGS_HPP
#define BAMR_FLAGS_HPP

#include "config.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
// Refinement status
////////////////////////////////////////////////////////////////////////////////
enum class bamr_status_t {
  none = 0,
  refine = 1,
  coarsen = -1
};

////////////////////////////////////////////////////////////////////////////////
// Refinement status
////////////////////////////////////////////////////////////////////////////////
struct bamr_flags_t {
  
  bamr_status_t flags[3] = {
    bamr_status_t::none,
    bamr_status_t::none,
    bamr_status_t::none
  };

  bamr_status_t operator[](int_t i) const { return flags[i]; }
  bamr_status_t & operator[](int_t i) { return flags[i]; }

  auto begin() const { return &flags[0]; }
  auto end() const { return &flags[3]; }


  void mark_refine_all() {
    flags[0] = bamr_status_t::refine;
    flags[1] = bamr_status_t::refine;
    flags[2] = bamr_status_t::refine;
  }
  
  void mark_coarsen_all() {
    flags[0] = bamr_status_t::coarsen;
    flags[1] = bamr_status_t::coarsen;
    flags[2] = bamr_status_t::coarsen;
  }

  bool has_refine() const {
    return
      flags[0] == bamr_status_t::refine || 
      flags[1] == bamr_status_t::refine ||
      flags[2] == bamr_status_t::refine;
  }

};



} // namespae


#endif
