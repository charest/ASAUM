#include "file_type.hpp"
#include "utils/string_utils.hpp"

namespace prl {

file_type_t get_file_type(const std::string & str )
{
  auto base = basename(str);
  auto i = base.rfind( '.', base.length() );
  size_t cnt = 0;
  while ( i != std::string::npos ) {
    auto ext = base.substr(i+1, base.length()-1);
    if ( ext == "g" || ext == "exo" || ext == "e-s" ) {
      if ( cnt != 0 ) return file_type_t::partitioned_exodus;
      else            return file_type_t::exodus;
    }
    base = base.substr(0, i);
    i = base.rfind( '.', base.length() );
    ++cnt;
  }
  return file_type_t::unknown;
}

}
