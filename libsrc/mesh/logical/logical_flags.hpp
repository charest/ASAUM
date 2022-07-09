#ifndef LOGICAL_FLAGS_HPP
#define LOGICAL_FLAGS_HPP

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// \brief The block boundary flags.
///
/// They are in a seperate struct since they are not always needed
////////////////////////////////////////////////////////////////////////////////
struct logical_flags_t 
{ 
  bool flags[6] = {false, false, false, false, false, false};
  
  bool is_boundary(bool is_hi, int_t dir) const
  { return flags[2*dir + is_hi]; }
};


} // prl


#endif
