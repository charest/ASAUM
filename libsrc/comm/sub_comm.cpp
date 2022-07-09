#include "sub_comm.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Constructor
////////////////////////////////////////////////////////////////////////////////
sub_comm_t::sub_comm_t(mpi_comm_t & comm, bool include) 
{
  int_t world_size = comm.size();
  
  // determine the ranks with data
  included_.clear();
  included_.resize(world_size);
  comm.all_gather(include, included_);

  // renumber the ranks
  rank_map_.clear();
  rank_map_.resize(world_size, -1);
  owns_ = false;
  for  (int_t r=0, i=0; r<world_size; ++r) {
    if (included_[r]) {
      rank_map_[r] = i;
      i++;
    }
    else {
      owns_ = true;
    }
  }

  // build a new sub-comm
  if (owns_)
    comm_ = new mpi_comm_t(comm, included_);
  else
    comm_ = &comm;

}

} // namespace
