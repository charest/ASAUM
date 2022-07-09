#include "config.hpp"
#include "ctm.hpp"

#include <cmath>

namespace prl {


////////////////////////////////////////////////////////////////////////////////
/// Helper transformation function
////////////////////////////////////////////////////////////////////////////////
void ctm_rotate(int_t val, int_t & sign, int_t & index)
{
  sign = (val < 0) ? -1 : 1;
  index = std::abs(val) - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Rotate a transformation vector
////////////////////////////////////////////////////////////////////////////////
void ctm_rotate(int_t * in, const int_t * rot, int_t ndims)
{
  int_t temp[3];
  int_t rot_sign, rot_index;

  for (int_t d=0; d<ndims; ++d) {
    ctm_rotate(rot[d], rot_sign, rot_index);
    temp[d] = rot_sign * in[rot_index];
  }

  for (int_t d=0; d<ndims; ++d)
    in[d] = temp[d];
}


////////////////////////////////////////////////////////////////////////////////
/// Constructors
////////////////////////////////////////////////////////////////////////////////
ctm_t::ctm_t(const int_t * ctm, int_t ndims, bool include_offsets) :
  ctm_t(ndims)
{
  if (include_offsets)
    std::copy_n(ctm, 2*ndims, data_.begin());
  else
    std::copy_n(ctm, ndims, data_.begin());
}

ctm_t::ctm_t(const ctm_t & ctm, bool include_offsets) :
  ctm_t(ctm.num_dims())
{
  if (include_offsets)
    std::copy(ctm.data_.begin(), ctm.data_.end(), data_.begin());
  else
    std::copy_n(ctm.data_.begin(), ndims_, data_.begin());
}
  
////////////////////////////////////////////////////////////////////////////////
/// Transform coordinates
////////////////////////////////////////////////////////////////////////////////
void ctm_t::transform(int_t & i, int_t & j, int_t & k) const
{
  int_t vals[] = {i, j, k};
  int_t s, index;

  ctm_rotate(data_[0], s, index);
  i = s*vals[index] + data_[3+index];

  ctm_rotate(data_[1], s, index);
  j = s*vals[index] + data_[3+index];
  
  ctm_rotate(data_[2], s, index);
  k = s*vals[index] + data_[3+index];
}

void ctm_t::transform(int_t & i, int_t & j) const
{
  int_t vals[] = {i, j};
  int_t s, index;

  ctm_rotate(data_[0], s, index);
  i = s*vals[index] + data_[2+index];

  ctm_rotate(data_[1], s, index);
  j = s*vals[index] + data_[2+index];
}

void ctm_t::transform(int_t & i) const
{
  auto s = (data_[0] < 0) ? -1 : 1;
  i = s*i + data_[1];
}
  
void ctm_t::rotate(bool * vals) const
{
  bool tmp[3];
  int_t s, index;

  for (int_t i=0; i<ndims_; ++i) {
    ctm_rotate(data_[i], s, index);
    tmp[i] = s && vals[index];
  }

  for (int_t i=0; i<ndims_; ++i)
    vals[i] = tmp[i];
}
  
////////////////////////////////////////////////////////////////////////////////
/// Is another transformation alignned
////////////////////////////////////////////////////////////////////////////////
bool ctm_t::is_aligned(const ctm_t b) const
{ 
  for (int_t i=0; i<ndims_; ++i)
    if (data_[i] != b.data_[i])
      return false;
  return true;
}
  
////////////////////////////////////////////////////////////////////////////////
/// Zero offsets
////////////////////////////////////////////////////////////////////////////////
void ctm_t::zero_offsets() {
  for (int_t i=0; i<ndims_; ++i)
    offset(i) = 0;
}
  

} // prl
