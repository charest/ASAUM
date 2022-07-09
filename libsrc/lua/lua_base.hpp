#ifndef LUA_BASE_HPP
#define LUA_BASE_HPP

// use lua
extern "C" {
  #include <lua.h>
  #include <lauxlib.h>
}

#include <iostream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

namespace prl {

//! \brief Use a shared pointer to the lua state.
//! Multiple objects may use the same lua state, so we don't want to 
//! destroy it unless all objects are destroyed.
using lua_state_ptr_t = std::shared_ptr<lua_State>;

////////////////////////////////////////////////////////////////////////////////
/// \brief A base class for several of the implemented objects.
/// This class mainly contains a lua state pointer which all derived
/// classes will use.  It also has some utility member functions.
////////////////////////////////////////////////////////////////////////////////
class lua_base_t {

protected:

  /// \brief The state pointer.
  lua_state_ptr_t state_;

public:

  /// \brief Default constructor.
  lua_base_t()
    : state_( luaL_newstate(), [](lua_State * s) { lua_close(s); } )
  {}
  
  /// \brief Copy constructor.
  lua_base_t(const lua_state_ptr_t & state) 
    : state_(state)
  {}
  
  /// \brief Return the raw state pointer.
  /// \remark Non-const version.
  auto state() { return state_.get(); }
  /// \brief Return the raw state pointer.
  /// \remark Const version.
  auto state() const { return state_.get(); }

  /// \brief Get the ith row of the stack.
  /// \param [in] i The row of the stack.
  /// \return A string with the type and value at the ith row.
  std::string get_row(int i) const;
 
  /// \brief Dump the stack to an output stream/
  /// \param [in,out] os The stream to output to.
  /// \return A reference to the output stream. 
  std::ostream& dump_stack(std::ostream& os) const;
  
  /// \brief Print the last row in the stack.
  void print_last_row() const;

  /// \brief the output operator.
  /// \param [in,out] os  The output stream.
  /// \param [in] s  The object whose stack to dump.
  /// \return The output stream.
  friend std::ostream& operator<<(std::ostream& os, const lua_base_t & s)  
  { 
    return s.dump_stack(os); 
  }
};

} // namespace

#endif
