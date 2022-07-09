#ifndef SIDE_SET_HPP
#define SIDE_SET_HPP

#include "config.hpp"

#include <map>
#include <string>
#include <vector>

namespace prl {

class mpi_comm_t;

struct side_set_info_t {
  int_t id;
  std::string label;

  void pack(std::vector<unsigned char> & buffer) const;
  void unpack(unsigned char const * & buffer);

  static void broadcast(
      std::map<int_t, side_set_info_t> & side_sets,
      mpi_comm_t & comm);

};


} // namespace

#endif
