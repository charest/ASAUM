/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some macros for raising errors.
////////////////////////////////////////////////////////////////////////////////
#ifndef ERRORS_HPP
#define ERRORS_HPP

#include "config.hpp"
#include "formatter.hpp"

// system includes
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace prl {

///////////////////////////////////////////////////////////////////
//! \brief Raised when a general runtime error occurs.
///////////////////////////////////////////////////////////////////
class exception_t : public std::runtime_error
{
  std::string loc_str_;
public:
  exception_t(const std::string & str, const std::string & loc) 
    : std::runtime_error(str), loc_str_(loc) {}

  std::string loc_str() const { return loc_str_; }
};

} // namespace

//! if exceptions are enabled, set the throw macro
//#define THROW_EXCEPTION(e) do { throw e; } while (0)
#define THROW_EXCEPTION(e) do { \
  std::cout << std::endl; \
  std::cout << e.loc_str() << std::endl; \
  std::cout << e.what() << std::endl; \
  std::cout << std::endl; \
  std::abort(); \
} while (0)

#define PRL_HERE \
  __FILE__ << " :: " << __FUNCTION__ << " :: L" << __LINE__


////////////////////////////////////////////////////////////////////////////////
//! \brief Raise a runtime error.
//! \param[in] the message to display
////////////////////////////////////////////////////////////////////////////////
#define THROW_ERROR(msg)                                                       \
  do {                                                                         \
    THROW_EXCEPTION(prl::exception_t(prl::formatter_t() << msg,                \
      prl::formatter_t() << PRL_HERE));                                        \
  } while(0)


////////////////////////////////////////////////////////////////////////////////
//! \brief Assert that something is true.
//! \param[in] cond  the condition to assert evaluates to true
//! \param[in] msg  the message to display
////////////////////////////////////////////////////////////////////////////////
#define assert_true(cond, msg)                                          \
  if ( ! (cond) )                                                       \
    THROW_ERROR("Assertion falied: " << msg)

#endif
