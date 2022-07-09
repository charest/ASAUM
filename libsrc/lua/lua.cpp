#include "lua.hpp"
#include "lua_modules.hpp"

namespace prl {

/// searcher for required modules
static int searcher(lua_State* const L)
{
  // Get the module name.
  const char* const modname = lua_tostring(L, 1);

  // search for the module
  auto it = Modules.find(modname);
      
  // Found the module
  if (it != Modules.end()) {
    const int res = luaL_loadfilex(L, it->second.c_str(), "t");
    if (res != LUA_OK) return lua_error(L);
    return 1;
  } // found

  // Oops...
  lua_pushfstring(L, "unknown module \"%s\"", modname);
  return 1;
}



/// \brief Main constructor.
/// \param [in] with_system  If true, load all system libraries.  
///                          Default is true.
lua_t::lua_t(bool with_system, bool with_extras) : lua_base_t() 
{
  if (!state_)
    THROW_ERROR("Cannot initialize lua state.");
  // open all system libraries
  if ( with_system )
    luaL_openlibs(state());

  if (with_extras) register_searcher();
}

/// \brief Run a string through the lua interpreter.
/// \param [in] script  The script to run.
/// \return The lua error code.
bool lua_t::run_string( const std::string & script )
{
  auto ret = luaL_dostring(state(),script.c_str());
  if ( ret ) {
    print_last_row();
    THROW_ERROR("Cannot load buffer.");
  }
  return ret;
}

/// \brief Load a file in the lua interpreter.
/// \param [in] file  The file to load.
void lua_t::loadfile( const std::string & file )
{
  auto ret = luaL_dofile(state(),file.c_str());
  if ( ret ) {
    print_last_row();
    THROW_ERROR("Cannot load file.");
  }
}

/// \brief Access an object in the global table.
/// \param [in] key  The key to access.
/// \return A lua_result_t object which points to the value of the table 
///         lookup.
lua_result_t lua_t::operator[]( const std::string & key ) const &
{
  auto s = state();
  // the function name
  lua_getglobal(s, key.c_str());
  // get the size of the object
  auto len = lua_rawlen(s, -1);
  // return the global object with a pointer to a location in the stack
  return { state_, key, make_lua_ref(state_), len };
}

/// register the extras search function
void lua_t::register_searcher() {
  auto s = state();
  // Get the package global table.
  lua_getglobal(s, "package");
  // Get the list of searchers in the package table.
  lua_getfield(s, -1, "searchers");
  // Get the number of existing searchers in the table.
  const size_t length = lua_rawlen(s, -1);

  // Add our own searcher to the list.
  lua_pushcfunction(s, searcher);
  lua_rawseti(s, -2, length + 1);

  // Remove the seachers and the package tables from the stack.
  lua_pop(s, 2);
}


} // namespace
