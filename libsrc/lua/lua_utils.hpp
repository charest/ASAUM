#ifndef LUA_UTILS_HPP
#define LUA_UTILS_HPP

// user includes
#include "lua.hpp"
#include "lua_result.hpp"
#include "utils/string_utils.hpp"
#include "utils/result.hpp"

// system libraries
#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <functional>
#include <memory>
#include <sstream>
#include <tuple>
#include <vector>

namespace prl {


#define lua_try_access_as(state, key, ...)                                     \
  (!state[key].empty()) ?                                                      \
    state[key].template as<__VA_ARGS__>() :                                             \
    throw std::runtime_error(                                                  \
      "\'" + state[key].name() +                                               \
      "\' does not exist in the lua state you are accessing."                  \
    )

#define lua_try_access(state, key)                                             \
  (!state[key].empty()) ?                                                      \
    state[key] :                                                               \
    throw std::runtime_error(                                                  \
      "\'" + state[key].name() +                                               \
      "\' does not exist in the lua state you are accessing."                  \
    )


////////////////////////////////////////////////////////////////////////////////
/// Does the input exist?
////////////////////////////////////////////////////////////////////////////////
bool validate(lua_t res, const std::string & key, bool is_verbose);

////////////////////////////////////////////////////////////////////////////////
/// Does the input exist?
////////////////////////////////////////////////////////////////////////////////
bool validate(lua_result_t res, const std::string & key, bool is_verbose);

////////////////////////////////////////////////////////////////////////////////
/// Does the input exist and match the provided options?
////////////////////////////////////////////////////////////////////////////////
template<typename T>
bool validate(
    lua_result_t res,
    const std::string & key,
    const std::vector<T> & options,
    bool is_verbose)
{
  if (!validate(res, key, is_verbose)) return false;
  auto val = res[key].as<T>();
  for (const auto & opt : options)
    if (opt == val) return true;
  if (is_verbose) {
    std::cout << "ERROR: Unknown value of '" << val << "' for key '"
      << key << "' under '" << res.name() << "'." << std::endl
      << "Acceptable options include: ";
    for (const auto & opt : options)
      std::cout << opt << " ";
    std::cout << std::endl;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Does the input exist and is it representable as the correct type
////////////////////////////////////////////////////////////////////////////////
template<typename T>
std::vector<T> as_vector(
    lua_result_t top,
    const std::string & key,
    bool is_verbose)
{

  using vector_t = std::vector<T>;

  if (!validate(top, key, is_verbose)) return {};

  auto res = top[key];

  auto is_table = res.is_table();

  auto ty = typeid(T).name();

  // vector
  if (is_table) {
    if (res.is_representable_as<vector_t>())
      return res.as<vector_t>();
    else {
      if (is_verbose) {
        std::cout << "ERROR: Cannot represent values for "  << top.name()
          << "["  << key << "] as {" << ty << ", " << ty << ", ...}"
          << std::endl;
      }
      return {};
    }
  }
  // scalar
  else {
    if (res.is_representable_as<T>())
      return {res.as<T>()};
    else {
      if (is_verbose) {
        std::cout << "ERROR: Cannot represent values for "  << top.name()
          << "["  << key << "] as " << ty
          << std::endl;
      }
      return {};
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Does the input exist and is it representable as the correct type and size
////////////////////////////////////////////////////////////////////////////////
template<typename T>
std::vector<T> as_vector(
    lua_result_t top,
    const std::string & key,
    bool is_verbose,
    std::size_t size)
{

  using vector_t = std::vector<T>;

  if (!validate(top, key, is_verbose)) return {};

  auto res = top[key];

  auto is_table = res.is_table();
  
  auto ty = typeid(T).name();

  // vector
  if (is_table) {
    if (res.size() != size) {
      if (is_verbose) {
        std::cout << "ERROR: Expected " << size << " values for " << top.name()
          << "["  << key << "], found " << res.size() << std::endl;
      }
      return {};
    }
    if (res.is_representable_as<vector_t>())
      return res.as<vector_t>();
    else {
      if (is_verbose) {
        std::cout << "ERROR: Cannot represent values for "  << top.name()
          << "["  << key << "] as {" << ty << ", " << ty << ", ...}"
          << std::endl;
      }
      return {};
    }
  }
  // scalar
  else {
    if (size != 1) {
      if (is_verbose) {
        std::cout << "ERROR: Expected " << size << " values for " << top.name()
          << "["  << key << "], found 1" << std::endl;
      }
      return {};
    }
    if (res.is_representable_as<T>())
      return {res.as<T>()};
    else {
      if (is_verbose) {
        std::cout << "ERROR: Cannot represent values for "  << top.name()
          << "["  << key << "] as " << ty
          << std::endl;
      }
      return {};
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Does the input exist and is it representable as the correct type and size
////////////////////////////////////////////////////////////////////////////////
template<typename T>
result_t<T> as_scalar(
    lua_result_t top,
    const std::string & key,
    bool is_verbose)
{

  if (!validate(top, key, is_verbose)) return {false};

  auto res = top[key];

  auto is_table = res.is_table();
  
  auto ty = typeid(T).name();

  // vector
  if (is_table) {
    if (is_verbose) {
      std::cout << "ERROR: Expected scalar value for " << top.name()
        << "["  << key << "], found vector of length " << res.size() << std::endl;
    }
    return {false};
  }
  // scalar
  else {
    if (res.is_representable_as<T>())
      return {res.as<T>(), true};
    else {
      if (is_verbose) {
        std::cout << "ERROR: Cannot represent value for "  << top.name()
          << "["  << key << "] as " << ty
          << std::endl;
      }
      return {false};
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Does the input exist and is it representable as the correct type and sizie
////////////////////////////////////////////////////////////////////////////////
bool validate_function(
    lua_result_t top,
    const std::string & key,
    bool is_verbose);

} // namespace prl
  
#endif
