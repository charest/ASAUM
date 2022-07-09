#ifndef CAST_HPP
#define CAST_HPP

#include <algorithm>

namespace prl {

//----------------------------------------------------------------------------
// A utilitiy for casting to byte arrays
//----------------------------------------------------------------------------
template<typename DATA_TYPE, typename BUFFER_TYPE>
void
cast_insert(DATA_TYPE const * const data, size_t len, BUFFER_TYPE & buf) {
  using buf_value_t = std::decay_t<decltype(buf[0])>;
  auto n = sizeof(DATA_TYPE) * len;
  auto p = reinterpret_cast<const buf_value_t *>(data);
  buf.insert(buf.end(), p, p + n);
}

//----------------------------------------------------------------------------
// A utilitiy for casting to byte arrays
//----------------------------------------------------------------------------
template<typename DATA_TYPE, typename BUFFER_TYPE>
void
cast(DATA_TYPE const * const data, size_t len, BUFFER_TYPE *& buf) {
  auto n = sizeof(DATA_TYPE) * len;
  auto p = reinterpret_cast<BUFFER_TYPE const * const>(data);
  std::copy(p, p + n, buf);
  buf += n;
}

//----------------------------------------------------------------------------
// A utilitiy for casting from byte arrays
//----------------------------------------------------------------------------
template<typename BUFFER_TYPE, typename DATA_TYPE>
void
uncast(BUFFER_TYPE const *& buffer, size_t len, DATA_TYPE * data) {
  auto n = sizeof(DATA_TYPE) * len;
  using data_type_no_cv = typename std::remove_cv<DATA_TYPE>::type;
  auto data_no_cv = const_cast<data_type_no_cv *>(data); // because of std::string
  auto p = reinterpret_cast<BUFFER_TYPE *>(data_no_cv);
  std::copy(buffer, buffer + n, p);
  buffer += n;
}

} // namespace

#endif
