#ifndef RESULT_HPP
#define RESULT_HPP

namespace prl {

template<typename T>
struct result_t {
  T value;
  bool exists;
  result_t(bool ex) : value(), exists(ex) {}
  result_t(const T & val, bool ex) : value(val), exists(ex) {}
  operator bool() const { return exists; }
};


}

#endif
