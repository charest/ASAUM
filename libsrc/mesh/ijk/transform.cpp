#include "config.hpp"
#include "ctm.hpp"
#include "transform.hpp"

namespace prl {


////////////////////////////////////////////////////////////////////////////////
/// Constructors
////////////////////////////////////////////////////////////////////////////////

transform_t::transform_t(int_t ndims) :
  ndims_(ndims), tm_(ndims*ndims, 0), offset_(ndims, 0)
{
  for (int_t i=0; i<ndims; ++i)
    tm_[i*ndims + i] = 1;
}

transform_t::transform_t(const int_t * ctm, int_t ndims, bool include_offsets)
  : ndims_(ndims), tm_(ndims*ndims, 0), offset_(ndims, 0)
{
  for (int_t i=0; i<ndims; ++i)
    this->operator()( std::abs(ctm[i])-1, i) = sign(ctm[i]);
  
  if (include_offsets)
    std::copy_n(ctm+ndims, ndims, offset_.begin());
}
  
transform_t::transform_t(const ctm_t & ctm) :
  transform_t(&ctm.data_[0], ctm.num_dims(), true)
{}

transform_t::transform_t(const transform_t & t1, const transform_t & t2) :
  transform_t( std::min(t1.num_dims(), t2.num_dims()) )
{
  for (int_t i=0; i<ndims_; ++i) {
    for (int_t j=0; j<ndims_; ++j) {
      tm_[i*ndims_+j]=0;
      auto & entry = tm_[i*ndims_+j];
      entry = 0;
      for (int_t k=0; k<ndims_; ++k)
        entry += t1[i*ndims_+k]*t2[k*ndims_+j];
    }
  }
}
  
////////////////////////////////////////////////////////////////////////////////
/// Transose a matrix
////////////////////////////////////////////////////////////////////////////////
void transform_t::transpose() {
  for (int_t i=0; i<ndims_; ++i)
    for (int_t j=i+1; j<ndims_; ++j)
      std::swap(tm_[i*ndims_+j], tm_[j*ndims_+i]);
}
  
////////////////////////////////////////////////////////////////////////////////
/// Transform coordinates
////////////////////////////////////////////////////////////////////////////////
void transform_t::transform(int_t & ii, int_t & jj, int_t & kk) {
  int_t i=ii, j=jj, k=kk;
  ii = tm_[0]*i + tm_[1]*j + tm_[2]*k + offset_[0];
  jj = tm_[3]*i + tm_[4]*j + tm_[5]*k + offset_[1];
  kk = tm_[6]*i + tm_[7]*j + tm_[8]*k + offset_[2];
}

void transform_t::transform(int_t & ii, int_t & jj) {
  int_t i=ii, j=jj;
  ii = tm_[0]*i + tm_[1]*j + offset_[0];
  jj = tm_[2]*i + tm_[3]*j + offset_[1];
}

void transform_t::transform(int_t & ii) 
{ ii = tm_[0]*ii + offset_[0]; }


////////////////////////////////////////////////////////////////////////////////
/// Reverse a transformation
////////////////////////////////////////////////////////////////////////////////
void transform_t::reverse() {
  
  int_t offset[3];
  for (int_t i=0; i<ndims_; ++i)
    offset[i] = -offset_[i];
  
  transpose();
  std::memset(offset_.data(), 0, offset_.size()*sizeof(int_t));
  transform(offset);

  for (int_t i=0; i<ndims_; ++i)
    offset_[i] = offset[i];
}

////////////////////////////////////////////////////////////////////////////////
/// Create a compact transformation from this matrix
////////////////////////////////////////////////////////////////////////////////
void transform_t::compact(int_t * ctm, bool include_offsets) const
{
  for (int_t i=0; i<ndims_; ++i) {
    ctm[i] = 0;
    for (int_t j=0; j<ndims_; ++j)
      ctm[i] += (j+1)*tm_[j*ndims_+i];
  }

  if (include_offsets) {
    for (int_t i=0; i<ndims_; ++i)
      ctm[ndims_+i] = offset_[i];
  }
}
  
void transform_t::compact(ctm_t & ctm) const
{ compact(&ctm.data_[0], true); }

////////////////////////////////////////////////////////////////////////////////
/// Compact a transformation and return it
////////////////////////////////////////////////////////////////////////////////
ctm_t compact(const transform_t & t)
{
  ctm_t ctm(t.num_dims());
  t.compact(ctm);
  return ctm;
}


////////////////////////////////////////////////////////////////////////////////
/// Transpose a transformation and return it
////////////////////////////////////////////////////////////////////////////////
transform_t transpose(const transform_t & t)
{
  auto tmp = t;
  tmp.transpose();
  return tmp;
}

////////////////////////////////////////////////////////////////////////////////
/// Reverse a transformation and return it
////////////////////////////////////////////////////////////////////////////////
transform_t reverse(const transform_t & t)
{
  auto tmp = t;
  tmp.reverse();
  return tmp;
}
  
////////////////////////////////////////////////////////////////////////////////
/// Output operator
////////////////////////////////////////////////////////////////////////////////
std::ostream &operator<<( std::ostream &out, const transform_t & trans )
{
  ctm_t ctm(trans.num_dims());
  trans.compact(ctm);
  out << ctm;
  return out;
}

} // prl
