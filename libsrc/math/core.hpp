#ifndef CORE_HPP
#define CORE_HPP

namespace prl {

//! 3d dot product
inline real_t dot_3d(const real_t * a, const real_t * b)
{ 
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

//! 2d dot product
inline real_t dot_2d(const real_t * a, const real_t * b)
{ 
  return a[0]*b[0] + a[1]*b[1];
}

//! nd dot product
inline real_t dot_nd(const real_t * a, const real_t * b, int_t n)
{ 
  real_t dot=0;
  for (int_t i=0; i<n; ++i)
    dot += a[i]*b[i];
  return dot;
}

//! sgn function
inline int sgn(real_t val) {
  return (0 < val) - (val < 0);
}

//! compute the square of a number
inline real_t sqr(real_t x) { return x*x; }

}

#endif
