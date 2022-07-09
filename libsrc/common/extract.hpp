#ifndef EXTRACT_HPP
#define EXTRACT_HPP

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// extract the  entries out of an array
////////////////////////////////////////////////////////////////////////////////
template<
  typename T,
  typename U,
  typename V
>
void extract(
    const T & arr,
    const U & ids,
    std::vector<V> & res,
    int_t ncomp = 1)
{
  auto num = ids.size();
  res.clear();
  res.resize(ncomp*num);
  for (size_t i=0; i<num; ++i) {
    auto id = ids[i];
    for (int_t j=0; j<ncomp; ++j)
      res[i*ncomp + j] = arr[id*ncomp + j];
  }
}

} // namesapce

#endif
