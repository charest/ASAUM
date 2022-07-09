#include "arguments.hpp"
#include "errors.hpp"
#include "string_utils.hpp"

#include <algorithm>

namespace prl {

namespace detail {
constexpr arg_info_t type_tr<bool>::value;
constexpr arg_info_t type_tr<float>::value;
constexpr arg_info_t type_tr<double>::value;
constexpr arg_info_t type_tr<std::string>::value;
}
  
///////////////////////////////////////////////////////////////////////////////
// output operators
///////////////////////////////////////////////////////////////////////////////
std::ostream & operator<<(std::ostream & out, const option_description_t & obj)
{
  out << obj.desc_;
  return out;
}

std::ostream & operator<<(std::ostream & out, const value_description_t & obj)
{
  out << obj.desc_;
  return out;
}

///////////////////////////////////////////////////////////////////////////////
//! \brief Option base intializer
///////////////////////////////////////////////////////////////////////////////
void option_base_t::init() {
  auto & args = arguments_t::instance();
  if (name_.empty())
    args.add_positional(this);
  else
    args.add_option(this);
}
  
///////////////////////////////////////////////////////////////////////////////
/// Is this a valid option argument?
///////////////////////////////////////////////////////////////////////////////
bool option_base_t::is_valid() const {
  if (info_.isbool) return true;
  if (results_.empty()) return false;
  for (const auto & i : results_)
    if (!(*info_.validate)(i))
      return false;
  return true;
}

///////////////////////////////////////////////////////////////////////////////
/// Add an option
///////////////////////////////////////////////////////////////////////////////
void arguments_t::add_option(option_base_t * opt)
{
  if (opt->name_.empty())
    std::cout << "Warning: Empty option." << std::endl;
  
  // add the category
  auto res = category_set_.emplace(opt->cat_);
  if (res.second)
    categories_.emplace_back(opt->cat_);
  
  category_options_[opt->cat_].emplace_back(opt);
    
  options_.emplace_back(opt);
  
  auto tks = split(opt->name_, {' ', ',', ';'});
  for (const auto & tk : tks) {
    auto res = option_map_.emplace(tk, opt);
    if (!res.second) {
      std::cout << "Command line option with token '" << tk;
      std::cout << "' already exists!" << std::endl;
      THROW_ERROR("Duplicate command-line option!");
    }
  }

}

///////////////////////////////////////////////////////////////////////////////
/// Add a positional argument
///////////////////////////////////////////////////////////////////////////////
void arguments_t::add_positional(option_base_t * opt)
{ positional_.emplace_back(opt); }

///////////////////////////////////////////////////////////////////////////////
/// Find a matchnig argument
///////////////////////////////////////////////////////////////////////////////
option_base_t * arguments_t::find(const std::string & arg) 
{ 
  auto it = option_map_.find(arg);
  if (it != option_map_.end())
    return it->second;
  return nullptr;
}
  
///////////////////////////////////////////////////////////////////////////////
/// Insert a dash
///////////////////////////////////////////////////////////////////////////////
std::string arguments_t::dash(const std::string & str)
{
  if (str.length() == 1)
    return "-";
  else
    return "--";
}
  
///////////////////////////////////////////////////////////////////////////////
/// get the value description string
///////////////////////////////////////////////////////////////////////////////
std::string arguments_t::value_string(option_base_t * opt)
{
  std::stringstream ss;
  std::string val_desc;
  
  if (opt->value_desc_)
    val_desc = opt->value_desc_.desc_;
  else if (!opt->info_.isbool)
    val_desc = to_uppercase(opt->info_.str);

  if (!val_desc.empty()) {
    if (opt->info_.islist) {
      ss << "<" << val_desc << ",..." << ">";
    }
    else {
      ss << "<" << val_desc << ">";
    }
  }

  return ss.str();
}
  
///////////////////////////////////////////////////////////////////////////////
/// Print the category usage
///////////////////////////////////////////////////////////////////////////////
void arguments_t::print_category(const option_category_t * cat)
{

  //----------------------------------------------------------------------------
  // First get max width
  
  std::vector< std::pair<std::string, const option_base_t*> > single_strms;
  std::vector< std::pair<std::string, const option_base_t*> > double_strms;

  size_t maxwidth = 0;
  for (auto ptr : category_options_.at(cat)) {
    
    std::stringstream ss;
    
    // split the name
    auto tks = split(ptr->name_, {' ', ',', ';'});
    const auto & tk = tks[0];
    
    ss << " " << dash(tk) << tk;

    if (tks.size()>1) {
      ss << " [";
      auto & tk = tks[1];
      ss << dash(tk) << tk;
      for (size_t i=2; i<tks.size(); ++i) {
        auto & tk = tks[i];
        ss << ", " << dash(tk) << tk;
      }
      ss << "]";
    }
  
    ss << " " << value_string(ptr);

    maxwidth = std::max( maxwidth, ss.str().length() );

    if (tk.length() == 1)
      single_strms.emplace_back( std::make_pair(ss.str(), ptr) );
    else
      double_strms.emplace_back( std::make_pair(ss.str(), ptr) );

  }
  
  //----------------------------------------------------------------------------
  // Now print
  for (const auto & strm : single_strms) {
    std::cout << std::setw(maxwidth) << std::left << strm.first;
    std::cout << " - " << strm.second->desc_ << std::endl;
  }

  //if (double_strms.size()) std::cout << std::endl;
  
  for (const auto & strm : double_strms) {
    std::cout << std::setw(maxwidth) << std::left << strm.first;
    std::cout << " - " << strm.second->desc_ << std::endl;
  }

  std::cout << std::endl;
}
  
  
///////////////////////////////////////////////////////////////////////////////
/// Print the usage
///////////////////////////////////////////////////////////////////////////////
void arguments_t::print_usage() {

  std::cout << "USAGE: " << appname_;
  for (auto ptr : options_) {
    if (ptr->required_) {
      auto tks = split(ptr->name_, {' ', ',', ';'});
      const auto & tk = tks[0];
      std::cout << " " << dash(tk) << tk;
      std::cout << " " << value_string(ptr);
    }
  }
  std::cout << " [options]";
  for (auto opt : positional_) {
    auto req = opt->required_;
    std::cout << " ";
    if (!req) std::cout << "[";
    std::cout << value_string(opt);
    if (!req) std::cout << "]";
  }
  std::cout << std::endl << std::endl;
  
  size_t maxwidth = 0;
  std::stringstream pos_ss;
  for (auto opt : positional_) {
    pos_ss << " " << value_string(opt);
    maxwidth = std::max( maxwidth, pos_ss.str().length() );
  }
  
  for (auto opt : positional_) {
    std::cout << std::setw(maxwidth) << std::left << value_string(opt);
    std::cout << " - " << opt->desc_ << std::endl;
  }
  
  if (positional_.size()) std::cout << std::endl;

  
  for (auto cat : categories_) {
    if (cat) 
      std::cout << cat->desc_ << std::endl; 
    else
      std::cout << "OPTIONS: " << std::endl;

    print_category(cat);
  }

}

///////////////////////////////////////////////////////////////////////////////
/// Main argument parsing
///////////////////////////////////////////////////////////////////////////////
int arguments_t::parse(int argc, char ** argv, bool verbose)
{

  // store args
  appname_ = argv[0];
  for (int i=1; i<argc; ++i)
    args_.emplace_back(argv[i]);

  std::size_t positional_counter = 0;
  bool force_positional = false;
  option_base_t * current = nullptr;
  auto nargs = args_.size();

  //----------------------------------------------------------------------------
  // now loop through args and process them
  for (size_t i=0; i<nargs; ++i) {

    // search for arg
    auto & arg = args_[i];
  
    if (!force_positional && arg.size()==2 && arg=="--")
    {
      force_positional = true;
      current = nullptr;
      continue;
    }

    auto isarg = arg.size()>1 && arg[0] == '-' &&
      !std::isdigit(arg[1]) && !force_positional;

    //----------------------------------
    // arg was found, is it good?
    if (isarg) {
      
      // the search string
      int pos = 1;
      if (arg.size()>1 && arg[1] == '-')
        pos++;
      auto search = arg.substr(pos); 

      // search for it
      auto it = find(search);
      if (it) {
        if (search == "h" || search == "help") {
          if (verbose) print_usage();
          return 1;
        }
        else {
          it->present_ = true;
          if (it->info_.isbool) {
            it->results_.emplace_back("1");
            current = nullptr;
          }
          else {
            current = it;
          }
        }
      }
      else {
        if (verbose) {
          std::cout << "Unknown option '" << arg << "'." << std::endl;
          std::cout << std::endl;
          print_usage();
        }
        return -1;
      }
    
    }
    //----------------------------------
    // An arg value
    else if (current) {
      auto & res = current->results_;
      auto islist = current->info_.islist;
      if (!islist && res.size()) {
        if (verbose) {
          std::cout << "Option '" << dash(current->name_);
          std::cout << current->name_ << "' is not a list,";
          std::cout << " at least 2 detected.";
          std::cout << std::endl << std::endl;
          print_usage();
        }
        return -1;
      }
      else {
        res.emplace_back(arg);
        if (!islist) current = nullptr;
      }
    }
    //----------------------------------
    // Positional
    else if (positional_.size()) {
      if (positional_counter >= positional_.size()) {
        if (verbose) {
          std::cout << "Too many positional arguments.  ";
          std::cout << "Declared " << positional_.size();
          std::cout << " but at least " << positional_counter+1;
          std::cout << " detected.";
          std::cout << std::endl << std::endl;
          print_usage();
        }
        return -1;
      }

      auto pos = positional_[positional_counter];
      auto & res = pos->results_;
      auto islist = pos->info_.islist;

      if (!islist && res.size()) {
        if (verbose) {
          std::cout << "Positional argument is not a list, ";
          std::cout << "at least 2 detected.";
          std::cout << std::endl << std::endl;
          print_usage();
        }
        return -1;
      }
      else {
        pos->present_ = true;
        res.emplace_back(arg);
        if (pos->info_.isbool) positional_counter++;
      }
    }
    //----------------------------------
    // Error, no positional
    else {
      if (verbose) {
        std::cout << "Argument '" << arg << "' is positional, but no";
        std::cout << " positional arguments were declared.";
        std::cout << std::endl << std::endl;
        print_usage();
      }
      return -1;
    } // arg


  } // nargs
  
  //----------------------------------------------------------------------------
  // make sure required args are present and valid
  for (const auto ptr : options_) {
    
    // check presence
    if (ptr->required_ && !ptr->present_) 
    {
      if (verbose) {
        std::cout << "Required argument '" << dash(ptr->name_);
        std::cout << ptr->name_ << "' " << "is not present.";
        std::cout << std::endl << std::endl;
        print_usage();
      }
      return -1;
    }

    // check validity
    if (ptr->present_ && !ptr->is_valid()) {
      if (verbose) {
        std::cout << "Arguments for '" << dash(ptr->name_);
        std::cout << ptr->name_;
        std::cout << "' are not present or representable ";
        std::cout << "as " << ptr->info_.str << "s";
        std::cout << std::endl << std::endl;
        print_usage();
      }
      return -1;
    }
  
  } // nargs
  
  //----------------------------------------------------------------------------
  // make sure required positional are present and valid
  for (auto pos : positional_) {

    // check presence
    if (pos->required_ && !pos->present_) 
    {
      if (verbose) {
        std::cout << "Required positional argument is not present.";
        std::cout << std::endl << std::endl;
        print_usage();
      }
      return -1;
    }

    // check validity
    if (pos->present_ && !pos->is_valid()) {
      if (verbose) {
        std::cout << "Positional arguments are present or not representable ";
        std::cout << "as " << pos->info_.str << "s.";
        std::cout << std::endl << std::endl;
        print_usage();
      }
      return -1;
    }

  } // positional
    
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
//! \brief Parse the command line options
///////////////////////////////////////////////////////////////////////////////
int parse_command_line_arguments(int argc, char ** argv, bool verbose)
{
  auto & args = arguments_t::instance();
  option_t<bool> help(
      "h,help",
      option_description_t("Display avaialble options"));
  return args.parse(argc, argv, verbose);
}

} // nammespace
