#ifndef VECTOR_UTILS_HPP
#define VECTOR_UTILS_HPP

#include <vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Remove a list of ids from a vector
/// \note ids must be sorted
////////////////////////////////////////////////////////////////////////////////
template<typename T, typename U>
void remove(std::vector<T> & vec, std::vector<U> & ids)
{

  if (ids.empty()) return;

  auto it = ids.begin();

  size_t i = 0;
  for (size_t j=0; it != ids.end(); ++j) {

    // dont erase
    if (j != static_cast<size_t>(*it)) {
      vec[i] = vec[j]; 
      ++i;
    }
    // erase
    else {
      ++it;
    }

  } // for

  vec.resize(i);

}


} // namespace


#endif
