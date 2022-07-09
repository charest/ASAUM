#ifndef LUA_REF_HPP
#define LUA_REF_HPP

#include "lua_base.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// \brief A class to keep track of a lua reference.
////////////////////////////////////////////////////////////////////////////////
class lua_ref_t : public lua_base_t {

  /// \brief The reference is stored as a shared_ptr to an int.
  /// A shared pointer is used so that the reference isn't deleted until
  /// all associated objects are destroyed.  Multiple copies of a single
  /// reference may exist.
  std::shared_ptr<int> ref_;

  /// \brief Also store the type id of the object we are referencing
  //int type_ = LUA_TNONE;

public:

  /// Delete the default destructor.
  lua_ref_t() = delete;

  /// \brief The main constructor.
  /// References are created in LUA_REGISTRYINDEX table.
  /// \param [in] state  A pointer to a lua state.
  /// \param [in] ref  The lua reference key.
  lua_ref_t ( const lua_state_ptr_t & state, int ref, int /*type*/ )
    : lua_base_t(state), 
      ref_( 
        new int{ref},
        [s=state](int * r) 
        {
          luaL_unref(s.get(), LUA_REGISTRYINDEX, *r);
          delete r;
        } 
      )
      //type_(type)
  {}

  /// \brief Constructor to create an empty reference.
  lua_ref_t( const lua_state_ptr_t & state )
    : lua_ref_t(state, LUA_REFNIL, LUA_TNONE)
  {}

  /// \brief Push the refered value onto the stack.
  void push() const
  {
    lua_rawgeti(state(), LUA_REGISTRYINDEX, *ref_);
  }

  /// \brief return true if the pointed reference is null
  bool empty() const
  {
    return (*ref_ == LUA_REFNIL || *ref_ == LUA_NOREF);
  }

};

/// \brief Create a lua reference to the last value on the stack.
/// \param [in] state A lua state pointer.
/// \return A new lua_ref_t object.
inline lua_ref_t make_lua_ref(const lua_state_ptr_t & state)
{
  auto s = state.get();
  return { state, luaL_ref(s, LUA_REGISTRYINDEX), lua_type(s, -1) };
}


} // namepsace

#endif
