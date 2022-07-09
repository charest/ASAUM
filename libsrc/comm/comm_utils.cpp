#include "mpi_comm.hpp"

#include <string>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// exchange a string
/// \remark this is only needed if you don't know the size of the vector
////////////////////////////////////////////////////////////////////////////////
mpi_request_t broadcast_string(
    mpi_comm_t & comm,
    std::string & str,
    int src)
{
  size_t size;
  bool is_root = comm.rank() == src;
  if (is_root) size = str.size();
  comm.broadcast(size, src);
  if (!is_root) {
    str.clear();
    str.resize(size);
  }
  return comm.broadcast(str, src);
}

} // namespace
