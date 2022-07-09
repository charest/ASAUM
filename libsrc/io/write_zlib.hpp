#ifndef WRITE_ZLIB_HPP
#define WRITE_ZLIB_HPP

#include "config.hpp"

#include <vector>

namespace prl {

std::vector<byte_t> zlib_deflate(const byte_t * bytes_to_encode, std::size_t in_len);


} // namespace

#endif
