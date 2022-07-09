////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Utilities for string operations.
////////////////////////////////////////////////////////////////////////////////
#ifndef STRING_UTILS_HPP
#define STRING_UTILS_HPP

// system includes
#include <cstring>
#include <iomanip>
#include <iterator>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

namespace prl {


//! \brief replace all occuraces of "from" to "to"
std::string replace_all(std::string str, const std::string & from, const std::string & to);

////////////////////////////////////////////////////////////////////////////////
//! \brief Convert a value to a string.
//! \param [in] x The value to convert to a string.
//! \return the new string
////////////////////////////////////////////////////////////////////////////////
template < typename T >
auto to_string(const T & x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

template < typename T >
auto to_string(const std::vector<T> & x) {
  std::stringstream ss;
  ss << "(";
  for (std::size_t i=0; i<x.size()-1; ++i)
    ss << x[i] << ", ";
  ss << x.back() << ")";
  return ss.str();
}


//! \brief Get a file name.
std::string basename(const std::string & str); 

//! \brief Get the extension of a file name.
std::string file_extension(const std::string & str); 

//! \brief Remove the extension from a filename
std::string remove_extension(const std::string & str);

//! \brief Remove the extension from a filename
std::pair<std::string, std::string>
split_extension(const std::string & str);

//! \brief split a string using a list of delimeters
std::vector<std::string> split(
  const std::string & str, 
  std::vector<char> delim = {' '});

//! \brief the number of digits in a number
unsigned num_digits( unsigned n );

//! \brief Tack on an iteration number to a string
std::string zero_padded(unsigned int n, unsigned int padding = 6);

//! \brief reset a string stream
void reset(std::stringstream & ss);

//! check if its a number
bool is_integer(const std::string &);
bool is_signed_integer(const std::string &);
bool is_float(const std::string &);

//! convert a string to upper case
std::string to_uppercase(std::string str);

////////////////////////////////////////////////////////////////////////////////
/// hash a string
////////////////////////////////////////////////////////////////////////////////
template<typename T, typename U>
inline constexpr T
string_hash(U && str, const std::size_t n) {
  if(n == 0)
    return 0;

  // String-to-integer hash function, based on prime numbers.
  // References:
  //    https://stackoverflow.com/questions/8317508/hash-function-for-a-string
  //    https://planetmath.org/goodhashtableprimes

  const T P = 3145739; // prime
  const T Q = 6291469; // prime, a bit less than 2x the first

  T h = 37; // prime
  for(std::size_t i = 0; i < n; ++i)
    h = (h * P) ^ (str[i] * Q);
  return h;
} // string_hash


} // namespace

#endif
