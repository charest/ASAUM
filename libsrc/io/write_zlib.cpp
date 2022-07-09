#include "write_zlib.hpp"

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#include <iostream>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// compress with zlib
////////////////////////////////////////////////////////////////////////////////
std::vector<byte_t> zlib_deflate(const byte_t* bytes_to_encode, std::size_t in_len) {

  std::vector<byte_t> res;

#ifdef HAVE_ZLIB

  constexpr int CHUNK = 0x4000;

  z_stream strm;
  strm.zalloc = Z_NULL;
  strm.zfree  = Z_NULL;
  strm.opaque = Z_NULL;

  // see https://zlib.net/manual.html
  auto ret = deflateInit(
      &strm,
      Z_DEFAULT_COMPRESSION);
  if (ret < 0) {
    std::cout << "deflateInit2 returned a bad status of " << ret << std::endl;
    return {};
  }


  strm.next_in = (unsigned char *)bytes_to_encode;
  strm.avail_in = in_len;

  unsigned char out[CHUNK];

  do {
    int have;
    strm.avail_out = CHUNK;
    strm.next_out = out;
    auto ret = deflate(&strm, Z_FINISH);
    if (ret < 0) {
      std::cout << "deflate returned a bad status of " << ret << std::endl;
      return {};
    }
    have = CHUNK - strm.avail_out;
    res.insert(res.end(), out, out+have);

  } while (strm.avail_out == 0);
  
  deflateEnd (& strm);


#endif

  return res;

}

} // namespace
