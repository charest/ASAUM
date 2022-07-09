#ifndef CTM_HPP
#define CTM_HPP

#include "config.hpp"
#include "range.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

namespace prl {

class transform_t;
  
////////////////////////////////////////////////////////////////////////////////
/// Helper transformation function
////////////////////////////////////////////////////////////////////////////////
void ctm_rotate(int_t val, int_t & sign, int_t & index);

////////////////////////////////////////////////////////////////////////////////
/// Rotate a transformation vector
////////////////////////////////////////////////////////////////////////////////
void ctm_rotate(int_t * in, const int_t * rot, int_t ndims);
  

////////////////////////////////////////////////////////////////////////////////
/// Compact transformation
////////////////////////////////////////////////////////////////////////////////
class ctm_t {

  int_t ndims_ = 0;
  std::vector<int_t> data_;

public:

  ctm_t(int_t ndims) : ndims_(ndims), data_(2*ndims, 0)
  { std::iota(&data_[0], &data_[ndims], 1); }
  
  ctm_t(const int_t * ctm, int_t ndims, bool include_offsets=false);
  
  ctm_t(const ctm_t & ctm, bool include_offsets);

  int_t num_dims() const { return ndims_; }

  int_t & operator[](int_t i) { return data_[i]; }
  int_t operator[](int_t i) const { return data_[i]; }
  
  int_t & offset(int_t i) { return data_[ndims_+i]; }
  int_t offset(int_t i) const { return data_[ndims_+i]; }

  int_t * data() { return &data_[0]; }
  const int_t * data() const { return &data_[0]; }

  void zero_offsets();
  
  void rotate(bool * vals) const;
  
  template<typename T>
  void rotate(T * vals) const
  {
    T tmp[3];
    int_t s, index;

    for (int_t i=0; i<ndims_; ++i) {
      ctm_rotate(data_[i], s, index);
      tmp[i] = s*vals[index];
    }

    for (int_t i=0; i<ndims_; ++i)
      vals[i] = tmp[i];
  }
  
  template<typename T>
  void rotate_abs(T * vals) const
  {
    T tmp[3];
    int_t s, index;

    for (int_t i=0; i<ndims_; ++i) {
      ctm_rotate(data_[i], s, index);
      tmp[i] = vals[index];
    }

    for (int_t i=0; i<ndims_; ++i)
      vals[i] = tmp[i];
  }
  
  template<typename T>
  void scale1(T * levels)
  {
    for (int_t i=0; i<ndims_; ++i) {
      offset(i) *= std::pow(2, levels[i]);
      if (data_[i] < 0) {
        for (int_t l=0, f=1; l<levels[i]; ++l, f*=2)
          offset(i) += f;
      }
    }
  }
  
  template<typename T>
  void scale2(T * levels)
  {
    for (int_t i=0; i<ndims_; ++i)
      offset(i) *= std::pow(2, levels[i]);
  }
  

  template<typename T>
  void transform(T * vals) const
  {
    T tmp[3];
    int_t s, index;

    for (int_t i=0; i<ndims_; ++i) {
      ctm_rotate(data_[i], s, index);
      tmp[i] = s*vals[index] + data_[ndims_+index];
    }

    for (int_t i=0; i<ndims_; ++i)
      vals[i] = tmp[i];
  }
  
  void transform(int_t & i, int_t & j, int_t & k) const;
  
  void transform(int_t & i, int_t & j) const;
  
  void transform(int_t & i) const;

  template<typename T>
  void transform(T * begin, T * end) const
  { 
    for (int_t i=0; i<ndims_; ++i)
      end[i]--;
    transform(begin);
    transform(end);

    for (int_t i=0; i<ndims_; ++i) {
      if (data_[i]<0) std::swap(begin[i], end[i]);
      end[i]++;
    }
  }

  void transform(range_ijk_t & range) const
  { transform(&range.begin[0], &range.end[0]); }

  bool operator==(const ctm_t & b) const
  { return data_ == b.data_; }

  bool operator!=(const ctm_t & b) const
  { return data_ != b.data_; }

  bool is_aligned(const ctm_t b) const;
  
  friend std::ostream &operator<<( std::ostream &out, const ctm_t & ctm )
  {
    auto ndims = ctm.num_dims();
    
    out << "or=(";
    for (int_t i=0; i<ndims-1; ++i)
      out << ctm[i] << ", ";
    out << ctm[ndims-1] << ") ";

    out << "os=(";
    for (int_t i=0; i<ndims-1; ++i) 
      out << ctm.offset(i) << ", ";
    out << ctm.offset(ndims-1) << ")";

    return out;
  }

  friend class transform_t;

};


} // namespace

#endif
