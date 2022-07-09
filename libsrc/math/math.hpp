#ifndef MATH_HPP
#define MATH_HPP

#include <vector>

namespace prl {

//! compute the square of a number
template <typename T>
T sqr(T x) { return x*x; }

//! check the bounds
template<typename T>
bool bounded(T val, T minval, T maxval)
{ return val > minval && val < maxval; }

template<typename T>
bool bounded(const std::vector<T> val, T minval, T maxval)
{ 
  for (auto v : val) {
    if (v < minval || v > maxval) {
      return false;
    }
  }
  return true;
}

template <typename T>
int sign(T val)
{
  return (T(0) < val) - (val < T(0));
}

} // namespace

#endif
