#ifndef UTILS_TRACE_HPP
#define UTILS_TRACE_HPP 

#include "config.hpp"
#include "errors.hpp"

#include <chrono>
#include <iostream>
#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace prl {

using chrono_t = std::chrono::steady_clock;
using time_point_t = std::chrono::time_point<chrono_t>;
using duration_t = std::chrono::duration<double>;

class mpi_comm_t;

////////////////////////////////////////////////////////////////////////////////
/// Gathered tracer node linked list
////////////////////////////////////////////////////////////////////////////////
class gathered_trace_node_t {

  std::string name_;

  std::map<std::string, std::unique_ptr<gathered_trace_node_t>> children_;

  std::map<int_t, double> elapsed_;
  std::map<int_t, size_t> counter_;

public:

  void unpack(const byte_t *& begin, int_t rank);
  
  void print(size_t level, double tot) const;
  void dump(std::ostream & out, size_t level, double tot, size_t nranks) const;

  double elapsed() const;
  size_t counts() const;
  
  std::string pretty_name(size_t len) const;
};


////////////////////////////////////////////////////////////////////////////////
/// Tracer node linked list
////////////////////////////////////////////////////////////////////////////////
class trace_node_t {

  std::string name_ = "__top__";

  trace_node_t * parent_ = nullptr;
  std::unordered_map<std::string, std::unique_ptr<trace_node_t>> children_;

  duration_t elapsed_ = std::chrono::seconds(0);
  time_point_t start_;
  size_t counter_ = 0;

public:

  trace_node_t() { start(); }
  
  trace_node_t(const std::string & name, trace_node_t * parent)
    : name_(name), parent_(parent)
  {}

  trace_node_t * child(const std::string & name)
  { 
    auto it = children_.find(name);
    if (it == children_.end()) {
      auto child = std::make_unique<trace_node_t>(name, this);
      auto res = children_.emplace(name, std::move(child));
      return res.first->second.get();
    }
    else {
      return it->second.get();
    }
  }

  const std::string & name() const { return name_; }
  std::string pretty_name(size_t len) const;

  trace_node_t * parent() const { return parent_; }

  double elapsed() const { return elapsed_.count(); }

  void start() {
    start_ = chrono_t::now();
    counter_++;
  }

  void stop() { elapsed_ += chrono_t::now() - start_; }

  void check(const std::string & name)
  {
    if (name != name_)
      THROW_ERROR(
          "Trace annotation mismatch.  Expecting to be at '" << name <<
          "', but actually at '" << name_ << "'." );
  }

  void count(size_t & sz) const;
  void pack(std::vector<byte_t> & buff) const;
  
  void print(size_t level, double tot) const;
};

////////////////////////////////////////////////////////////////////////////////
/// Main tracer singleton
////////////////////////////////////////////////////////////////////////////////
class tracer_t {

  trace_node_t root_;
  trace_node_t * current_ = &root_;

public:

  static tracer_t & instance() {
    static tracer_t inst;
    return inst;
  }
  
  void start(const std::string & name)
  { 
    current_ = current_->child(name);
    current_->start();
  }

  void stop(const std::string & name)
  {
    current_->check(name);
    current_->stop();
    current_ = current_->parent();
  }

  void dump(mpi_comm_t & comm, bool print, const std::string & filename = "");

private:
  tracer_t() = default;
  tracer_t(const tracer_t&)= delete;
  tracer_t& operator=(const tracer_t&)= delete;
  
};

////////////////////////////////////////////////////////////////////////////////
/// Main trace class that whose scope is used to time a region
////////////////////////////////////////////////////////////////////////////////
class scoped_trace_t {

  std::string name_;

public:

  scoped_trace_t(const char * name) : name_(name) 
  { tracer_t::instance().start(name); }
  
  scoped_trace_t(const std::string & name) : scoped_trace_t(name.c_str()) {}

  ~scoped_trace_t() { tracer_t::instance().stop(name_); }

};

} // namespace

#endif
