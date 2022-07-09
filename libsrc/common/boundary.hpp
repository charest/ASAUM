#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "config.hpp"

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace prl {

class lua_result_t;
class mpi_comm_t;

////////////////////////////////////////////////////////////////////////////////
/// Boundary info
////////////////////////////////////////////////////////////////////////////////
class boundary_t {
  using function_type = std::function<bool(const std::vector<real_t> &)>;
  int_t id_ = -1;
  function_type where_;
  std::string name_;
public:
  boundary_t(const std::string & name) : name_(name) {}
  boundary_t(int_t id, const std::string & name) : id_(id), name_(name) {}
  const std::string & name() const { return name_; }
  int_t id() const { return id_; }
  void where(function_type fun) { where_ = fun; }
  function_type where() const { return where_; }
  bool has_where() const { return static_cast<bool>(where_); }
  virtual ~boundary_t() {}
};

} // namespace

#endif
