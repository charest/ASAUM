#ifndef UTILS_ARGUMENTS_HPP
#define UTILS_ARGUMENTS_HPP

#include "string_utils.hpp"
#include "type_traits.hpp"

// system includes
#include<iomanip>
#include<iostream>
#include<map>
#include<set>
#include<sstream>
#include<vector>


namespace prl {

class arguments_t;

//! \brief Parse the command line options
int parse_command_line_arguments(int argc, char ** argv, bool verbose);

namespace detail {

///////////////////////////////////////////////////////////////////////////////
// check if strings are valid
///////////////////////////////////////////////////////////////////////////////
inline bool converts(const std::string &)
{ return true; }

typedef bool (*validate_fun_t)(const std::string &);

///////////////////////////////////////////////////////////////////////////////
// Type string struct 
///////////////////////////////////////////////////////////////////////////////

/// the argument type info type
struct arg_info_t {
  const char * str;
  bool isbool;
  bool islist;
  validate_fun_t validate;
};

/// specializations to statically determine the argument info
template<typename T, typename U = void>
struct type_tr;

template<>
struct type_tr<bool>
{
  static constexpr arg_info_t value = arg_info_t{"bool", true, false, &converts};
  static bool convert(const std::string & str)
  { return !str.empty(); }
};

template<typename T>
struct type_tr<
  T,
  std::enable_if_t<
    std::is_integral<T>::value &&
    std::is_signed<T>::value &&
    !std::is_same<T, bool>::value
  >
>
{
  static constexpr arg_info_t value = arg_info_t{"int", false, false, &is_integer};
  static T convert(const std::string & str)
  { return std::stoi(str); }
};

template<typename T>
struct type_tr<
  T,
  std::enable_if_t<std::is_integral<T>::value && !std::is_signed<T>::value>
>
{
  static constexpr arg_info_t value = arg_info_t{"uint", false, false,
    &is_signed_integer};
  static T convert(const std::string & str)
  { return std::stoi(str); }
};

template<>
struct type_tr<float>
{
  static constexpr arg_info_t value = arg_info_t{"float", false, false,
    &is_float};
  static float convert(const std::string & str)
  { return std::stof(str); }
};

template<>
struct type_tr<double>
{
  static constexpr arg_info_t value = arg_info_t{"double", false, false,
    &is_float};
  static double convert(const std::string & str)
  { return std::stod(str); }
};


template<>
struct type_tr<std::string>
{
  static constexpr arg_info_t value = arg_info_t{"string", false, false, converts};
  static const std::string & convert(const std::string & str)
  { return str; }
};

template<typename T>
struct type_tr<std::vector<T>>
{
  static constexpr arg_info_t value = arg_info_t{
    type_tr<T>::value.str,
    false, 
    true,
    type_tr<T>::value.validate
  };
  
  static std::vector<T> convert(const std::vector<std::string> & strs)
  { 
    auto n = strs.size();
    std::vector<T> ret(n);
    for (size_t i=0; i<n; ++i)
      ret[i] = type_tr<T>::convert(strs[i]);
    return ret;
  }
};

template<typename T>
constexpr arg_info_t type_tr<std::vector<T>>::value;

} // namespace  

///////////////////////////////////////////////////////////////////////////////
//! \brief An option category 
///////////////////////////////////////////////////////////////////////////////
struct option_category_t {
  std::string desc_;
  option_category_t() = default;
  option_category_t(const std::string & desc) : desc_(desc) {}
};

///////////////////////////////////////////////////////////////////////////////
//! \brief An option description 
///////////////////////////////////////////////////////////////////////////////
struct option_description_t {
  std::string desc_;
  option_description_t() = default;
  option_description_t(const std::string & desc) : desc_(desc) {}
  friend std::ostream & operator<<(std::ostream &, const option_description_t &);
};

///////////////////////////////////////////////////////////////////////////////
//! \brief An option value description 
///////////////////////////////////////////////////////////////////////////////
struct value_description_t {
  std::string desc_;
  value_description_t() = default;
  value_description_t(const std::string & desc) : desc_(desc) {}
  operator bool() const { return !desc_.empty(); }
  friend std::ostream & operator<<(std::ostream &, const value_description_t &);
};

///////////////////////////////////////////////////////////////////////////////
//! \brief Category description 
///////////////////////////////////////////////////////////////////////////////
struct category_description_t {
  std::string desc_;
  category_description_t() = default;
  category_description_t(const std::string & desc) : desc_(desc) {}
};

///////////////////////////////////////////////////////////////////////////////
/// enum class for optional/required variables
///////////////////////////////////////////////////////////////////////////////
enum class required_t {
  no,
  yes
};

///////////////////////////////////////////////////////////////////////////////
//! \brief An option
///////////////////////////////////////////////////////////////////////////////
struct option_base_t {

  bool present_ = false;

  std::string name_;
  detail::arg_info_t info_;

  const option_category_t * cat_ = nullptr;
  option_description_t desc_;
  value_description_t value_desc_;
  bool required_ = false;

  std::vector<std::string> results_;

  option_base_t(
      const std::string & name,
      detail::arg_info_t info,
      required_t req) :
    name_(name), info_(info), required_(req==required_t::yes)
  { init(); }
  
  option_base_t(
      const std::string & name,
      detail::arg_info_t info,
      const option_description_t & desc,
      required_t req) :
    name_(name), info_(info), desc_(desc), required_(req==required_t::yes) 
  { init(); }
  
  option_base_t(
      const std::string & name,
      detail::arg_info_t info,
      const option_description_t & desc,
      const value_description_t & val,
      required_t req) :
    name_(name), info_(info), desc_(desc), value_desc_(val),
    required_(req==required_t::yes)
  { init(); }
  
  option_base_t(
      const std::string & name,
      detail::arg_info_t info,
      const option_category_t & cat,
      required_t req) :
    name_(name), info_(info), cat_(&cat), required_(req==required_t::yes)
  { init(); }
  
  option_base_t(
      const std::string & name,
      detail::arg_info_t info,
      const option_category_t & cat,
      const option_description_t & desc,
      required_t req) :
    name_(name), info_(info), cat_(&cat), desc_(desc),
    required_(req==required_t::yes)
  { init(); }
  
  option_base_t(
      const std::string & name,
      detail::arg_info_t info,
      const option_category_t & cat,
      const option_description_t & desc,
      const value_description_t & val,
      required_t req) :
    name_(name), info_(info), cat_(&cat), desc_(desc), value_desc_(val),
    required_(req==required_t::yes)
  { init(); }

  void init();

  explicit operator bool() const
  { return present_; } 

  bool is_valid() const;
};


///////////////////////////////////////////////////////////////////////////////
//! \brief An option
///////////////////////////////////////////////////////////////////////////////
template<typename T>
struct option_t : public option_base_t {

  // regular arg constructors
  option_t(
      const std::string & name,
      required_t req = required_t::no) :
    option_base_t(name, traits(), req)
  {}
  
  option_t(
      const std::string & name,
      const option_description_t & desc,
      required_t req = required_t::no) :
    option_base_t(name, traits(), desc, req)
  {}
  
  option_t(
      const std::string & name,
      const option_description_t & desc,
      const value_description_t & val,
      required_t req = required_t::no) :
    option_base_t(name, traits(), desc, val, req)
  {}
  
  option_t(
      const std::string & name,
      const option_category_t & cat,
      required_t req = required_t::no) :
    option_base_t(name, traits(), cat, req)
  {}
  
  option_t(
      const std::string & name,
      const option_category_t & cat,
      const option_description_t & desc,
      required_t req = required_t::no) :
    option_base_t(name, traits(), cat, desc, req)
  {}
  
  option_t(
      const std::string & name,
      const option_category_t & cat,
      const option_description_t & desc,
      const value_description_t & val,
      required_t req = required_t::no) :
    option_base_t(name, traits(), cat, desc, val, req)
  {}
  
  // positional arg constructors
  option_t(
      required_t req = required_t::no) :
    option_base_t("", traits(), req)
  {}
  
  option_t(
      const option_description_t & desc,
      required_t req = required_t::no) :
    option_base_t("", traits(), desc, req)
  {}
  
  option_t(
      const option_description_t & desc,
      const value_description_t & val,
      required_t req = required_t::no) :
    option_base_t("", traits(), desc, val, req)
  {}
  
  option_t() = delete;
  option_t(const option_t &) = delete;
  option_t & operator=(const option_t &) = delete;
  option_t && operator=(const option_t &&) = delete;
 
  ~option_t() = default;
  
  static constexpr detail::arg_info_t traits()
  { return detail::type_tr<T>::value; }

  static constexpr bool is_bool()
  { return std::is_same<T, bool>::value; }
  
  static constexpr bool is_vec()
  { return is_vector<T>::value; }

  using option_base_t::operator bool;

  T value() { return *this; }

  template< bool Enabled = !is_bool() && !is_vec(),
    typename = typename std::enable_if_t<Enabled> >
  operator T() const
  { 
    if (present_) 
      return detail::type_tr<T>::convert(results_[0]);
    else
      return {};
  }
  
  template< bool Enabled = !is_bool() && is_vec(),
    typename std::enable_if_t<Enabled>* = nullptr >
  operator T() const
  {
    if (present_)
      return detail::type_tr<T>::convert(results_);
    else
      return {};
  }
  

};

///////////////////////////////////////////////////////////////////////////////
//! \brief The main arguments class
///////////////////////////////////////////////////////////////////////////////
class arguments_t {

  std::string appname_;
  std::vector<std::string> args_;

  std::set<const option_category_t*> category_set_;
  std::vector<const option_category_t*> categories_;

  std::map<const option_category_t*, std::vector<option_base_t*>> category_options_;

  std::vector<option_base_t*> options_;
  std::map<std::string, option_base_t*> option_map_;

  std::vector<option_base_t*> positional_;
  
  arguments_t() = default;
  ~arguments_t() = default;
  arguments_t(const arguments_t &) = delete;
  arguments_t & operator=(const arguments_t &) = delete;

  void add_option(option_base_t * opt);
  void add_positional(option_base_t * opt);

  option_base_t * find(const std::string & arg);
  
  void print_category(const option_category_t * cat);

public:
  
  static arguments_t & instance()
  {
    static arguments_t instance;
    return instance;
  }

  static std::string dash(const std::string & str);

  static std::string value_string(option_base_t * opt);

  void print_usage();
  
  int parse(int argc, char ** argv, bool verbose);

  friend struct option_base_t;
};

} // namespace

#endif
