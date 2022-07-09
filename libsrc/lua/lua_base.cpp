#include "lua_base.hpp"

namespace prl {

/// \brief Get the ith row of the stack.
/// \param [in] i The row of the stack.
/// \return A string with the type and value at the ith row.
std::string lua_base_t::get_row(int i) const
{
  auto s = state();
  std::stringstream os;
  auto t = lua_type(s, i);
  switch (t) {
    case LUA_TSTRING:  // strings
      os << "string >> " << lua_tostring(s, i);
      break;
    case LUA_TBOOLEAN:  // booleans
      os << "bool   >> " << (lua_toboolean(s, i) ? "true" : "false");
      break;
    case LUA_TNUMBER:  // numbers
      os << "number >> " << lua_tonumber(s, i);
      break;
    default:  // other values
      os << "other  >> " << lua_typename(s, t);
      break;
  }
  return os.str();
}

/// \brief Dump the stack to an output stream/
/// \param [in,out] os The stream to output to.
/// \return A reference to the output stream. 
std::ostream& lua_base_t::dump_stack(std::ostream& os) const 
{
  auto s = state();
  auto top = lua_gettop(s);
  if ( top ) {
    os << "Row : Type   >> Value" << std::endl;
    for (int i = 1; i <= top; i++) {  /* repeat for each level */
      os << std::setw(3) << i << " : ";
      os << get_row(i);
      os << std::endl;  // put a separator
    }
  }
  else {
    os << "(stack empty)" << std::endl;
  }
  return os;
}

/// \brief Print the last row in the stack.
void lua_base_t::print_last_row() const
{
  auto s = state();
  auto top = lua_gettop(s);
  std::cout << "Row : Type   >> Value" << std::endl;
  std::cout << std::setw(3) << top << " : " << get_row(-1) << std::endl;
  //lua_pop(state(), 1);
}


} // namespace
