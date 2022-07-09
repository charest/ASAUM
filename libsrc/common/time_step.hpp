#ifndef TIME_STEP_HPP
#define TIME_STEP_HPP

#include "config.hpp"

namespace prl {

struct time_step_t {
  real_t val = consts::real_max;
  long_t pos = -1;
  int flag = 0;
  time_step_t() = default;
  time_step_t(real_t dt) : val(dt) {}
};
  
inline bool compare_lt(const time_step_t & a, const time_step_t & b) 
{ return a.val < b.val; }

} // namespace

#endif
