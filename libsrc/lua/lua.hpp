#ifndef LUA_HPP
#define LUA_HPP

#include "lua_base.hpp"
#include "lua_ref.hpp"
#include "lua_result.hpp"

// use lua
extern "C" {
  #include <lualib.h>
}

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// \brief The top level object for the lua interface.
/// This is the object the user will instantiate.
////////////////////////////////////////////////////////////////////////////////
class lua_t : public lua_base_t {

public:

  /// \brief Main constructor.
  /// \param [in] with_system  If true, load all system libraries.  
  ///                          Default is true.
  lua_t(bool with_system = true, bool with_extras = false);

  /// \brief Run a string through the lua interpreter.
  /// \param [in] script  The script to run.
  /// \return The lua error code.
  bool run_string( const std::string & script );

  /// \brief Load a file in the lua interpreter.
  /// \param [in] file  The file to load.
  void loadfile( const std::string & file );

  /// \brief Access an object in the global table.
  /// \param [in] key  The key to access.
  /// \return A lua_result_t object which points to the value of the table 
  ///         lookup.
  lua_result_t operator[]( const std::string & key ) const &;

  /// \brief Run a string through the lua interpreter.
  /// \param [in] script  The script to run.
  /// \return The lua error code.
  auto operator()( const std::string & script )
  {
    return run_string( script );
  }

  /// register the extras searcher
  void register_searcher();
};

} // namespace

#endif
