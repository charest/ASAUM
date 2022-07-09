#include "side_set.hpp"
#include "comm/mpi_comm.hpp"
#include "utils/cast.hpp"
#include "utils/errors.hpp"

#include <numeric>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Pack/unpack side info
////////////////////////////////////////////////////////////////////////////////
void side_set_info_t::pack(std::vector<byte_t> & buffer) const
{
  cast_insert( &id, 1, buffer );
  size_t len = label.size();
  cast_insert( &len, 1, buffer );
  cast_insert( label.c_str(), len, buffer );
}

void side_set_info_t::unpack(byte_t const * & buffer)
{
  uncast( buffer, 1, &id );
  size_t len;
  uncast( buffer, 1, &len );
  label.resize(len);
  uncast( buffer, len, label.data() );
}

////////////////////////////////////////////////////////////////////////////////
/// Broadcast side set info
////////////////////////////////////////////////////////////////////////////////
void side_set_info_t::broadcast(
    std::map<int_t, side_set_info_t> & side_sets,
    mpi_comm_t & comm)
{

  auto num_side_sets = side_sets.size();

  size_t comm_size = comm.size();
  
  // now make sure everyone has all side set meta data
  size_t tot_side_sets{0};
  comm.all_reduce(
    num_side_sets,
    tot_side_sets,
    redop_t::sum);

  if ( tot_side_sets > 0 ) {

    // pack buffer
    std::vector<byte_t> sendbuf;
    for ( auto ss : side_sets ) {
      size_t side_id = ss.first;
      cast_insert( &side_id, 1, sendbuf );
      ss.second.pack(sendbuf);
    }

    // exchange buffer sizes
    int buf_len = sendbuf.size();
    std::vector<int> recvcounts(comm_size);
    comm.all_gather(buf_len, recvcounts);

    // compute receive displacements
    std::vector<int> recvdispls(comm_size+1);
    recvdispls[0] = 0;
    std::partial_sum(recvcounts.begin(), recvcounts.end(), &recvdispls[1]);
    
    // exchange data
    std::vector<byte_t> recvbuf(recvdispls[comm_size]);
    comm.all_gatherv(sendbuf, recvbuf, recvcounts, recvdispls);

    // unpack
    for ( size_t r=0; r<comm_size; ++r ) {
  
      auto start = recvdispls[r];
      auto end = recvdispls[r+1];
      const auto * buffer = &recvbuf[ start ];
      const auto * buffer_end = &recvbuf[ end ];
  
      for ( auto i=start; i<end; ) {

        const auto * buffer_start = buffer;

        // the side id
        size_t side_id;
        uncast( buffer, 1, &side_id );

        // the side info
        side_set_info_t new_ss;
        new_ss.unpack( buffer );
        side_sets.emplace( side_id, std::move(new_ss) );

        // increment pointer
        i += std::distance( buffer_start, buffer );

      }
      // make sre we are at the right spot in the buffer
      if (buffer != buffer_end)
        THROW_ERROR("Exodus side set unpacking mismatch.");
    }
  }
} // broadcast

} // namespace
