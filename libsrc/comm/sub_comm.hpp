#ifndef SUB_COMM_HPP
#define SUB_COMM_HPP

#include "config.hpp"
#include "mpi_comm.hpp"

#include <vector>

namespace prl {

class sub_comm_t {

  mpi_comm_t * comm_ = nullptr;

  bool owns_ = false;
  
  std::vector<char> included_;

  std::vector<int_t> rank_map_;

public:

  sub_comm_t(mpi_comm_t & comm, bool include); 

  ~sub_comm_t()
  { 
    if (owns_ && comm_) delete comm_;
    comm_ = nullptr;
  }

  auto & comm() { return *comm_; }

  auto map_rank(int_t r) { return rank_map_[r]; }
  bool is_included(int_t r) { return included_[r]; }

};

} // namespace

#endif
