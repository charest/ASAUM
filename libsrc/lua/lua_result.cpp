#include "lua_result.hpp"

namespace prl {
  
/// \brief Check the stack to make sure it can be expanded.
/// \param [in] extra The desired size to grow the stack by. 
void lua_result_t::check_stack(int extra) const
{
  auto s = state();
  auto ret = lua_checkstack(s, extra);
  if ( !ret ) {
    std::ostringstream ss;
    ss << "Cannot grow stack " << extra << " slots operating on element \""
       << name_ << "\"." << std::endl << "Current stack size is " 
       << lua_gettop(s) << ".";
    THROW_ERROR(ss.str());
  }
} 

void lua_result_t::check_nil(const std::string & name) const
{
  if ( lua_isnil(state(), -1) ) {
    print_last_row();
    THROW_ERROR("\"" << name << "\" does not exist.");
  }
}

void lua_result_t::check_table(const std::string & name) const
{
  check_nil(name);
  if ( !lua_istable(state(), -1) ) {
    print_last_row();
    THROW_ERROR("\"" << name << "\" is not a table.");
  }
}

/// \brief Check that a result is valid.
/// \parm [in] name  The name of the object we are checking.
void lua_result_t::check_result(const std::string & name) const
{
  if (lua_isnil(state(), -1)) {
    print_last_row();
    THROW_ERROR("\"" << name << "\" returned nil.");
  }
}

/// \brief Check that a result is a function.
/// \parm [in] name  The name of the object we are checking.
void lua_result_t::check_function(const std::string & name) const
{
  check_nil(name);
  if ( !lua_isfunction(state(), -1) ) {
    print_last_row();
    THROW_ERROR("\"" << name << "\" is not a function.");
  }
}

/// \brief Push the last value onto the stack.
/// \remark The stack references are collected in reverse order.
void lua_result_t::push_last() const
{
  check_stack(1);
  refs_.back().push();
}

/// \brief Push all referred values onto the stack.
/// \remark The stack references are collected in reverse order.
void lua_result_t::push_all() const
{
  // make sure the stack can handle this
  check_stack(refs_.size());
  // push everything onto the stack in reverse order
  for ( auto && r : refs_ ) std::forward<decltype(r)>(r).push();
}
  
/// \brief Return true if all references are nil.
bool lua_result_t::empty() const
{
  for ( const auto & r : refs_ )
    if (!r.empty()) return false;
  return true;
}

  
bool lua_result_t::is_table() const {
  auto s = state();
  push_last();
  auto res = lua_istable(s, -1);
  lua_remove(s,-1);
  return res;
}

bool lua_result_t::is_function() const {
  auto s = state();
  push_last();
  auto res = lua_isfunction(s, -1);
  lua_remove(s,-1);
  return res;
}
  
/// \brief Get the lua value in the table given a key.
/// \param [in] key  The key to search for.
/// \return A new lua_result_t with a reference to the result of the table
///         lookup.
lua_result_t lua_result_t::operator[]( const std::string & key ) const &
{
  auto s = state();
  auto new_name = name_+".[\""+key+"\"]";
  // push the table onto the stack
  push_last();
  // make sure we are accessing a table
  check_table(name_);
  // push the key onto the stack
  check_stack(1);
  lua_pushlstring(s, key.c_str(), key.size());
  // now get the table value, the key gets pushed from the stack
  lua_rawget(s, -2);
  // get the size of the object
  auto len = lua_rawlen(s, -1);
  // get a reference to the result, and pop the table from the stack
  auto ref = make_lua_ref(state_);
  lua_pop(s, -1);
  // return the global object with a pointer to a location in the stack
  return { state_, new_name, std::move(ref), len };
}

/// \brief Get the lua value in the table given an index.
/// \param [in] n  The index to access.
/// \return A new lua_result_t with a reference to the result of the table
///         lookup.
lua_result_t lua_result_t::operator[]( int n ) const &
{
  auto s = state();
  auto new_name = name_+"["+std::to_string(n)+"]";
  // push the table onto the stack
  push_last();
  // make sure we are accessing a table
  check_table(name_);
  // now get the table value, the key gets pushed from the stack
  lua_rawgeti(s, -1, n);
  // get the size of the object
  auto len = lua_rawlen(s, -1);
  // get a reference to the result, and pop the table from the stack
  auto ref = make_lua_ref(state_);
  lua_pop(s, -1);
  // return the global object with a pointer to a location in the stack
  return { state_, new_name, std::move(ref), len };
}


} // namespace
