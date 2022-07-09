#ifndef LUA_PUSH_HPP
#define LUA_PUSH_HPP

#include "utils/errors.hpp"

// use lua
extern "C" {
  #include <lua.h>
}

#include <sstream>
#include <string>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// \defgroup lua_push lua_push
/// \brief Functions to push values onto the lua stack.
////////////////////////////////////////////////////////////////////////////////
/// \{

/// \brief Push an integer onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] i  The integer to push.
inline std::size_t lua_push(lua_State * s, int i)
{
  lua_pushinteger( s, i );
  return 1;
}

/// \brief Push a size_t onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] i  The size_t to push.
inline std::size_t lua_push(lua_State * s, size_t i)
{
  lua_pushinteger( s, i );
  return 1;
}

/// \brief Push a long long onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] i  The long long to push.
inline std::size_t lua_push(lua_State * s, long long i)
{
  lua_pushinteger( s, i );
  return 1;
}

/// \brief Push a float onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] x  The float to push.
inline std::size_t lua_push(lua_State * s, float x)
{
  lua_pushnumber( s, x );
  return 1;
}

/// \brief Push a double onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] x  The double to push.
inline std::size_t lua_push(lua_State * s, double x) 
{
  lua_pushnumber( s, x );
  return 1;
}

/// \brief Push a boolean onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] b  The boolean to push.
inline std::size_t lua_push(lua_State * s, bool b)
{
  lua_pushboolean( s, b );
  return 1;
}

/// \brief Push a character array onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] str  The character array to push.
inline std::size_t lua_push(lua_State * s, const char * str)
{
  lua_pushstring( s, str );
  return 1;
}

/// \brief Push a std::string onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] str  The string to push.
inline std::size_t lua_push(lua_State * s, const std::string & str) 
{
  lua_pushlstring( s, str.c_str(), str.size() );
  return 1;
}

/// \brief Push a std::array onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] str  The array to push.
template <
  typename T,
  std::size_t N,
  template <typename, std::size_t> class Array
>
inline std::size_t lua_push(lua_State * s, const Array<T, N> & arr)
{
  // Caller always checks for 1, we must check for the rest
  auto ret = lua_checkstack(s, N - 1);
  if (!ret) {
    std::ostringstream ss;
    ss << "Cannot grow stack " << (N - 1) << " slots operating on element \""
       << "\"." << std::endl
       << "Current stack size is " << lua_gettop(s) << ".";
    THROW_ERROR(ss.str());
  }
  // push each element of the array
  lua_newtable(s);
  for (std::size_t i = 0; i < N; ++i) {
    lua_push(s, arr[i]);
    lua_rawseti(s,-2,i+1);
  }
  return 1;
}

/// \brief Push a std::vector onto the stack.
/// \param [in] s  The lua state to push a value to.
/// \param [in] str  The array to push.
template <
  typename T,
  typename...Args,
  template <typename, typename...> class Vector
>
inline std::size_t lua_push(lua_State * s, const Vector<T, Args...> & vec)
{
  // Caller always checks for 1, we must check for the rest
  auto N = vec.size();
  auto ret = lua_checkstack(s, N - 1);
  if (!ret) {
    std::ostringstream ss;
    ss << "Cannot grow stack " << (N - 1) << " slots operating on element \""
       << "\"." << std::endl
       << "Current stack size is " << lua_gettop(s) << ".";
    THROW_ERROR(ss.str());
  }
  // push each element of the array
  lua_newtable(s);
  for (std::size_t i = 0; i < N; ++i) {
    lua_push(s, vec[i]);
    lua_rawseti(s,-2,i+1);
  }
  return 1;
}
/// \}

}

#endif
