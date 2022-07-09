#ifndef LUA_VALUE_HPP
#define LUA_VALUE_HPP

#include "utils/errors.hpp"

// use lua
extern "C" {
  #include <lua.h>
}

#include <string>
#include <type_traits>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// \defgroup lua_value lua_value
/// \brief A struct to extract and typecast values from the lua stack.
////////////////////////////////////////////////////////////////////////////////
/// \{

/// \brief The default implementation.
/// \tparam T  The type.
template < typename T, typename Enable = void>
struct lua_value {};

/// \brief The implementation for integral values.
template <typename T>
struct lua_value< T, std::enable_if_t< std::is_integral<T>::value > > 
{
  /// \brief Return the value in the lua stack.
  /// \param [in] s  The lua state to query.
  /// \param [in] index  The row to access.  Defaults to the value at the top
  ///                    of the stack.
  /// \return The requested value as a long long.
  static T get(lua_State * s, int index = -1) 
  {
    if ( !lua_isnumber(s,index) )
     THROW_ERROR( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to int."
    );
    auto i = lua_tointeger(s, index);
    lua_remove(s,index);
    return i;
  }
  
  static bool is_representable(lua_State * s, int index = -1) 
  { 
    auto flag = lua_isnumber(s,index);
    lua_remove(s,index);
    return flag;
  }
};


/// \brief The implementation for double.
template <typename T>
struct lua_value<T, std::enable_if_t< std::is_floating_point<T>::value > >
{
  /// \brief Return the value in the lua stack.
  /// \param [in] s  The lua state to query.
  /// \param [in] index  The row to access.  Defaults to the value at the top
  ///                    of the stack.
  /// \return The requested value as a double.
  static T get(lua_State * s, int index = -1) 
  {
    if ( !lua_isnumber(s,index) )
     THROW_ERROR( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to double."
    );
    auto x = lua_tonumber(s, index);
    lua_remove(s,index);
    return x;
  }
  
  static bool is_representable(lua_State * s, int index = -1) 
  { 
    auto flag = lua_isnumber(s,index);
    lua_remove(s,index);
    return flag;
  }
};

/// \brief The implementation for bool.
template <>
struct lua_value<bool> 
{
  /// \brief Return the value in the lua stack.
  /// \param [in] s  The lua state to query.
  /// \param [in] index  The row to access.  Defaults to the value at the top
  ///                    of the stack.
  /// \return The requested value as a boolean.
  static bool get(lua_State * s, int index = -1) 
  {
    if ( !lua_isboolean(s,index) )
     THROW_ERROR( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to bool."
    );
    auto b = lua_toboolean(s, index);
    lua_remove(s, index);
    return b;
  }
  
  static bool is_representable(lua_State * s, int index = -1) 
  { 
    auto flag = lua_isboolean(s,index);
    lua_remove(s,index);
    return flag;
  }
};

/// \brief The implementation for std::string.
template <>
struct lua_value<std::string> 
{
  /// \brief Return the value in the lua stack.
  /// \param [in] s  The lua state to query.
  /// \param [in] index  The row to access.  Defaults to the value at the top
  ///                    of the stack.
  /// \return The requested value as a string.
  static std::string get(lua_State * s, int index = -1)
  {
    if ( !lua_isstring(s, index) )
     THROW_ERROR( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to string."
    );
    auto str = lua_tostring(s, index);
    lua_remove(s, index);
    return str;
  }
  
  static bool is_representable(lua_State * s, int index = -1) 
  { 
    auto flag = lua_isstring(s,index);
    lua_remove(s,index);
    return flag;
  }
};

/// \brief The implementation for vectors.
template< 
  typename T, 
  typename Allocator,
  template<typename,typename> class Vector
>
struct lua_value< Vector<T,Allocator> > 
{

  /// \brief Return the value in the lua stack.
  /// \param [in] s  The lua state to query.
  /// \param [in] index  The row to access.  Defaults to the value at the top
  ///                    of the stack.
  /// \return The requested value as a vector.
  static Vector<T,Allocator> get(lua_State * s, int index = -1)
  {
    // make sure we are accessing a table
    if ( !lua_istable(s, index) )
     THROW_ERROR( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to vector."
    );
    // get the size of the table
    auto n = lua_rawlen(s, -1);
    // extract the results
    Vector<T,Allocator> res;
    res.reserve(n);
    for ( unsigned i=1; i<=n; ++i ) {
      lua_rawgeti(s, -1, i);  // push t[i]
      res.emplace_back( lua_value<T>::get(s) );
    }
    // remove it from the stack
    lua_remove(s, index);
    return res;
  }
  
  static bool is_representable(lua_State * s, int index = -1) 
  { 
    auto flag = lua_istable(s,index);
    if (flag) {
      // get the size of the table
      auto n = lua_rawlen(s, -1);
      // extract the results
      for ( unsigned i=1; i<=n; ++i ) {
        lua_rawgeti(s, -1, i);  // push t[i]
        if (!lua_value<T>::is_representable(s)) {
          flag = false;
          break;
        }
      }
    }
    lua_remove(s,index);
    return flag;
  }
};

/// \brief The implementation for vectors.
template< 
  typename T, 
  std::size_t N, 
  template <typename,std::size_t> class Array 
>
struct lua_value< Array<T,N> > 
{

  /// \brief Return the value in the lua stack.
  /// \param [in] s  The lua state to query.
  /// \param [in] index  The row to access.  Defaults to the value at the top
  ///                    of the stack.
  /// \return The requested value as a vector.
  static Array<T,N> get(lua_State * s, int index = -1)
  {
    // make sure we are accessing a table
    // make sure we are accessing a table
    if ( !lua_istable(s, index) )
     THROW_ERROR( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to vector."
    );
    // get the size of the table
    auto n = lua_rawlen(s, -1);
    if ( n != N )
      THROW_ERROR( 
        "Expecting array of size"<<N<<", stack array is size " << n
      );
    // extract the results
    Array<T,N> res;
    int min_n = std::min<int>(n,N);
    for ( int i=1; i<=min_n; ++i ) {
      lua_rawgeti(s, -1, i);  // push t[i]
      res[i-1] = lua_value<T>::get(s);
    }
    // remove it from the stack
    lua_remove(s, index);
    return res;
  }
  
  static bool is_representable(lua_State * s, int index = -1) 
  { 
    auto flag = lua_istable(s,index);
    if (flag) {
      // get the size of the table
      auto n = lua_rawlen(s, -1);
      // extract the results
      for ( decltype(n) i=1; i<=n; ++i ) {
        lua_rawgeti(s, -1, i);  // push t[i]
        if (!lua_value<T>::is_representable(s)) {
          flag = false;
          break;
        }
      }
    }
    lua_remove(s,index);
    return flag;
  }
};

/// \}
  
/// \brief The implementation for ristra::matrix.
//! This lua_value_t expects to be initialized using nested lua tables,
//! e.g. {{1,2}, {3,4}}. Remember, ristra multiarrays are row-major ordered!
template< 
  typename T, 
  std::size_t M,
  std::size_t N,
  template <typename,std::size_t,std::size_t> class Matrix
>
struct lua_value< Matrix<T,M,N> > 
{

  /// \brief Return the value in the lua stack.
  /// \param [in] s  The lua state to query.
  /// \param [in] index  The row to access.  Defaults to the value at the top
  ///                    of the stack.
  /// \return The requested value as a matrix
  static Matrix<T,M,N> get(lua_State * s, int index = -1)
  {
    // make sure we are accessing a table
    if ( !lua_istable(s, index) )
     THROW_ERROR( "Invalid conversion of type \"" <<
      lua_typename(s, lua_type(s, index)) << "\" to multiarray."
    );
    // get the size of the table
    auto n = lua_rawlen(s, -1);
    if ( n != M )
      THROW_ERROR( 
        "Expecting array of size"<<M<<", stack array is size " << n
      );
    // extract the results
    Matrix<T,M,N> res;
    for ( int i=1; i<=M; ++i ) {
      lua_rawgeti(s, -1, i);  // push t[i]
      for (int j=1; j<=N; ++j) {
        lua_rawgeti(s, -1, j); // push t[j]
        res(i-1,j-1) = lua_value<T>::get(s);
      }
      lua_remove(s, -1);
    }
    // remove it from the stack
    lua_remove(s, index);
    return res;
  }
};

} // namespace

#endif
