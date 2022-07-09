// system includes
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

namespace prl {


////////////////////////////////////////////////////////////////////////////////
//! \brief replace all occuraces of "from" to "to"
//! \param [in] str  the input string
//! \param [in] from the string to search for
//! \param [in] to   the string to replace "from" with
//! \return the new string
////////////////////////////////////////////////////////////////////////////////
std::string
replace_all(std::string str, const std::string & from, const std::string & to) {
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
  }
  return str;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Get a file name.
//! \param [in] str  the input string
//! \return the base of the file name
////////////////////////////////////////////////////////////////////////////////
std::string basename(const std::string & str) 
{
#ifdef _WIN32
  char sep = '\\';
#else
  char sep = '/';
#endif

  auto i = str.rfind( sep, str.length() );
  if ( i != std::string::npos )
    return str.substr(i+1, str.length()-1);
  else
    return str;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Get the extension of a file name.
//! \param [in] str  the input string
//! \return the extension
////////////////////////////////////////////////////////////////////////////////
std::string file_extension(const std::string & str) 
{
  auto base = basename(str);
  auto i = base.rfind( '.', base.length() );

  if ( i != std::string::npos )
    return base.substr(i+1, base.length()-1);
  else
    return "";
 
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Remove the extension from a filename
//! \param [in] str  the input string
//! \return the name without extension
////////////////////////////////////////////////////////////////////////////////
std::string remove_extension(const std::string & str) {
    auto lastdot = str.find_last_of(".");
    if (lastdot == std::string::npos) return str;
    return str.substr(0, lastdot); 
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Remove the extension from a filename
//! \param [in] str  the input string
//! \return the name without extension
////////////////////////////////////////////////////////////////////////////////
std::pair<std::string, std::string>
split_extension(const std::string & str) {
  auto lastdot = str.find_last_of(".");
  if (lastdot == std::string::npos)
    return std::make_pair( str, "");
  else
    return std::make_pair( str.substr(0, lastdot), str.substr(lastdot+1) ); 
}


////////////////////////////////////////////////////////////////////////////////
//! \brief split a string using a list of delimeters
//! \param [in] str  the input string
//! \param [in] delim  the list of delimeters
//! \return the list of split strings
////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> split(
  const std::string & str, 
  std::vector<char> delim)
{

  if (str.empty()) return {};

  auto it = std::find(delim.begin(), delim.end(), ' ');
  bool has_space = (it != delim.end());

  const std::string * tmp_str = &str;
  if (!has_space) {
    auto tmp = new std::string;
    *tmp = replace_all(str, " ", "___SPACE___");
    tmp_str = tmp;
  }

  struct tokens_t : std::ctype<char>
  {
    using ctype_base = std::ctype_base;
    using cctype = std::ctype<char>;
    using ccmask = cctype::mask;
    
    tokens_t(const std::vector<char> & delims) 
      : cctype(get_table(delims)) {}

    static ctype_base::mask const * get_table(
      const std::vector<char> & delims
    ) {
      static const ccmask * const_rc = cctype::classic_table();
      static ccmask rc[cctype::table_size];
      std::memcpy(rc, const_rc, cctype::table_size*sizeof(ccmask));
      for (const auto & d : delims) 
        rc[static_cast<int>(d)] = ctype_base::space;
      return &rc[0];
    }
  };

  std::stringstream ss(*tmp_str);
  ss.imbue(std::locale(std::locale(), new tokens_t(delim)));
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> vstrings(begin, end);

  if (!has_space) {
    for (auto & s : vstrings)
      s = replace_all(s, "___SPACE___", " ");
    delete tmp_str;
  }

  return vstrings;
}

///////////////////////////////////////////////////////////////////////////////
/// The number of digits
///////////////////////////////////////////////////////////////////////////////
unsigned num_digits( unsigned n )
{
  auto number = n / 10;
  unsigned int digits = 1;
  while (number) {
    number /= 10;
    digits++;
  }

  return digits;

}


///////////////////////////////////////////////////////////////////////////////
//! \brief Tack on an iteration number to a string
///////////////////////////////////////////////////////////////////////////////
std::string zero_padded(unsigned int n, unsigned int padding)
{
  std::stringstream ss;
  ss << std::setw( padding ) << std::setfill( '0' ) << n;
  return ss.str();
}

///////////////////////////////////////////////////////////////////////////////
//! \brief reset a string stream
///////////////////////////////////////////////////////////////////////////////
void reset(std::stringstream & ss)
{
  ss.str( std::string() );
  ss.clear();
}

///////////////////////////////////////////////////////////////////////////////
//! check if its a number
///////////////////////////////////////////////////////////////////////////////
bool is_integer(const std::string & str)
{
  char* p;
  strtol(str.c_str(), &p, 10);
  return !*p;
}

bool is_signed_integer(const std::string & str)
{
  if (str.length() && str[0]=='-') return false;
  char* p;
  strtol(str.c_str(), &p, 10);
  return !*p;
}

///////////////////////////////////////////////////////////////////////////////
//! check if its a number
///////////////////////////////////////////////////////////////////////////////
bool is_float(const std::string & str)
{
  char* p;
  strtof(str.c_str(), &p);
  return !*p;
}

///////////////////////////////////////////////////////////////////////////////
//! convert to uppercase
///////////////////////////////////////////////////////////////////////////////
std::string to_uppercase(std::string str)
{
  std::for_each(str.begin(), str.end(), [](auto & a){a = ::toupper(a);});
  return str;
}


} // namespace
