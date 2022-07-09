#ifndef TRANSFORM_HPP
#define TRANSFORM_HPP

#include "config.hpp"
#include "range.hpp"
#include "math/math.hpp"

#include <array>
#include <cstring>

namespace prl {

class ctm_t;
class transform_t;

/// Compact a transformation and return it
ctm_t compact(const transform_t & t);

/// Transpose a transformation and return it
transform_t transpose(const transform_t & t);

/// Reverse a transformation and return it
transform_t reverse(const transform_t & t);


////////////////////////////////////////////////////////////////////////////////
/// Full transformation matrix
////////////////////////////////////////////////////////////////////////////////
class transform_t {

  int_t ndims_ = 0;
  std::vector<int_t> tm_;
  std::vector<int_t> offset_;

public:

  transform_t(const transform_t &) = default;
  
  transform_t(int_t ndims); 
  
  transform_t(const int_t * ctm, int_t ndims, bool include_offsets = false);
  
  transform_t(const ctm_t & ctm);
  
  transform_t(const transform_t & t1, const transform_t & t2);
  
  int_t num_dims() const { return ndims_; }

  void zero_offsets()
  { std::fill(offset_.begin(), offset_.end(), 0); }

  void transpose();
  
  int_t col_sum(int_t i) const {
    int_t sum = 0;
      for (int_t j=0; j<ndims_; ++j)
        sum += tm_[i*ndims_+j];
    return sum;
  }
  
  template<typename T>
  void rotate(T * vals)
  { 
    T tmp[3];
    
    for (int_t i=0; i<ndims_; ++i) {
      tmp[i] = 0;
      for (int_t j=0; j<ndims_; ++j)
        tmp[i] += tm_[i*ndims_+j] * vals[j];
    }

    for (int_t i=0; i<ndims_; ++i)
      vals[i] = tmp[i];
  }
  

  template<typename T>
  void transform(T * vals)
  { 
    T tmp[3];
    
    for (int_t i=0; i<ndims_; ++i) {
      tmp[i] = offset_[i];
      for (int_t j=0; j<ndims_; ++j)
        tmp[i] += tm_[i*ndims_+j] * vals[j];
    }

    for (int_t i=0; i<ndims_; ++i)
      vals[i] = tmp[i];
  }
  
  void transform(int_t & ii, int_t & jj, int_t & kk);
  
  void transform(int_t & ii, int_t & jj);

  void transform(int_t & ii);
  
  template<typename T>
  void transform(T * begin, T * end)
  { 
    for (int_t i=0; i<ndims_; ++i)
      end[i]--;
    transform(begin);
    transform(end);

    for (int_t i=0; i<ndims_; ++i) {
      auto sum = col_sum(i);
      if (sum<0) std::swap(begin[i], end[i]);
      end[i]++;
    }
  }

  void transform(range_ijk_t & range)
  { transform(&range.begin[0], &range.end[0]); }
  
  template<typename T>
  void scale1(T * levels)
  {
    for (int_t i=0; i<ndims_; ++i) {
      offset(i) *= std::pow(2, levels[i]);
      auto sum = col_sum(i);
      if (sum < 0) {
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
  

  void reverse();

  void compact(int_t * ctm, bool include_offsets=false) const;
  
  void compact(ctm_t & ctm) const;


  int_t & operator()(int_t i, int_t j)
  { return tm_[i*ndims_ + j]; }
  
  int_t operator()(int_t i, int_t j) const
  { return tm_[i*ndims_ + j]; }

  const int_t * offsets() const { return offset_.data(); }
  int_t * offsets() { return offset_.data(); }

  int_t offset(int_t i) const { return offset_[i]; }
  int_t & offset(int_t i) { return offset_[i]; }

private:

  int_t & operator[](int_t i)
  { return tm_[i]; }
  
  int_t operator[](int_t i) const
  { return tm_[i]; }

public:
  
  friend std::ostream &operator<<( std::ostream &out, const transform_t & trans );

};

} // namespace

#endif
