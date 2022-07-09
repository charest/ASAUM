#ifndef AMR_FLAGS_HPP
#define AMR_FLAGS_HPP

#include "config.hpp"

#include <ostream>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
// Refinement status
////////////////////////////////////////////////////////////////////////////////
enum class amr_status_t {
  none = 0,
  refine = 1,
  coarsen = -1
};

////////////////////////////////////////////////////////////////////////////////
// Refinement status
////////////////////////////////////////////////////////////////////////////////
struct amr_flags_t {
  
  amr_status_t flags[3] = {
    amr_status_t::none,
    amr_status_t::none,
    amr_status_t::none
  };

  amr_status_t operator[](int_t i) const { return flags[i]; }
  amr_status_t & operator[](int_t i) { return flags[i]; }

  auto begin() const { return &flags[0]; }
  auto end() const { return &flags[3]; }

  bool has_change() const {
    return
      flags[0] != amr_status_t::none || 
      flags[1] != amr_status_t::none ||
      flags[2] != amr_status_t::none;
  }

  void mark_refine(int_t d)
  { flags[d] = amr_status_t::refine; }

  void mark_coarsen(int_t d)
  { flags[d] = amr_status_t::coarsen; }

  void mark_refine_all() {
    flags[0] = amr_status_t::refine;
    flags[1] = amr_status_t::refine;
    flags[2] = amr_status_t::refine;
  }
  
  void mark_coarsen_all() {
    flags[0] = amr_status_t::coarsen;
    flags[1] = amr_status_t::coarsen;
    flags[2] = amr_status_t::coarsen;
  }

  bool has_refine() const {
    return
      flags[0] == amr_status_t::refine || 
      flags[1] == amr_status_t::refine ||
      flags[2] == amr_status_t::refine;
  }
  
  bool has_refine(int_t dim) const
  { return flags[dim] == amr_status_t::refine; }
  
  bool has_coarsen() const {
    return
      flags[0] == amr_status_t::coarsen || 
      flags[1] == amr_status_t::coarsen ||
      flags[2] == amr_status_t::coarsen;
  }
  
  bool has_coarsen(int_t dim) const
  { return flags[dim] == amr_status_t::coarsen; }
  
  friend std::ostream &operator<<( std::ostream &out, const amr_flags_t & f )
  {
    out << "(";
    for (int_t i=0; i<2; ++i)
      out << static_cast<int>(f[i]) << ", ";
    out << static_cast<int>(f[2]) << ") ";
    return out;
  }

};



} // namespae


#endif
