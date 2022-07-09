#ifndef FILE_TYPE_HPP
#define FILE_TYPE_HPP

#include <string>

namespace prl {

enum class file_type_t {
  exodus,
  partitioned_exodus,
  unknown
};

// try to figure out what kind of file it is
file_type_t get_file_type(const std::string & str );

} // namespace

#endif
