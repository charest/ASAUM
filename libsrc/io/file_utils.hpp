#ifndef FILE_UTILS_HPP
#define FILE_UTILS_HPP

namespace prl {

int make_dir(const char * name);

bool file_exists(const char * name);
int move_file(const char * src, const char * dst);

} // namespace

#endif
