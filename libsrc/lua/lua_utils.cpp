#include "lua_utils.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Does the input exist?
////////////////////////////////////////////////////////////////////////////////
bool validate(lua_t res, const std::string & key, bool is_verbose)
{
  if (res[key].empty()) {
    if (is_verbose) {
      std::cout << "ERROR: No input section for '" << key << "' exists."
        << std::endl
        << "A key/value pair '" << key << " = {...}' must exist in the top "
        << "scope of your input file."
        << std::endl;
    }
    return false;
  }
  return true;
}


////////////////////////////////////////////////////////////////////////////////
/// Does the input exist?
////////////////////////////////////////////////////////////////////////////////
bool validate(lua_result_t res, const std::string & key, bool is_verbose)
{
  if (res[key].empty()) {
    if (is_verbose) {
      std::cout << "ERROR: No input section for '" << key << "' exists under '"
        << res.name() << "'." << std::endl
        << "A key/value pair '" << key << " = {...}' must exist in your input file."
        << std::endl;
    }
    return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Does the input exist and is it representable as the correct type and sizie
////////////////////////////////////////////////////////////////////////////////
bool validate_function(
    lua_result_t top,
    const std::string & key,
    bool is_verbose)
{

  if (!validate(top, key, is_verbose)) return false;

  auto res = top[key];

  auto is_func = res.is_function();
  if (!is_func) {
    if (is_verbose) {
      std::cout << "ERROR: Expected function for " << top.name()
        << "["  << key << "]." << std::endl;
    }
    return false;
  }

  return true;
}

} // namespace
