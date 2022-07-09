#ifndef COMM_UTILS_HPP
#define COMM_UTILS_HPP

#include "mpi_comm.hpp"
#include "utils/crs.hpp"

#include <numeric>
#include <vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// exchange a vector
/// \remark this is only needed if you don't know the size of the vector
////////////////////////////////////////////////////////////////////////////////
template <typename T>
mpi_request_t broadcast_vector(
    mpi_comm_t & comm,
    std::vector<T> & vec,
    int src)
{
  size_t size;
  bool is_root = comm.rank() == src;
  if (is_root) size = vec.size();
  comm.broadcast(size, src);
  if (!is_root) {
    vec.clear();
    vec.resize(size);
  }
  return comm.broadcast(vec, src);
}
  
/// exchange a string
mpi_request_t broadcast_string(
    mpi_comm_t & comm,
    std::string & str,
    int src);
  
////////////////////////////////////////////////////////////////////////////////
/// determine offsets by exchanging counts
////////////////////////////////////////////////////////////////////////////////
template <typename T>
void global_offsets(
    mpi_comm_t & comm,
    const T & num,
    std::vector<T> & displs)
{
  auto comm_size = comm.size();
  
  displs.clear();
  displs.resize(comm_size+1);
  
  auto counts = make_array_ref( &displs[1], comm_size );
  comm.all_gather(num, counts);
  
  displs[0] = 0;
  std::partial_sum(displs.begin(), displs.end(), displs.begin());
}

////////////////////////////////////////////////////////////////////////////////
/// determine offsets by exchanging counts
////////////////////////////////////////////////////////////////////////////////
template <
  typename SendType,
  typename IdType,
  bool Enabled = is_container_v<SendType> &&
    is_container_v<IdType> &&
    std::is_same<int, typename IdType::value_type>::value, 
  typename = typename std::enable_if_t<Enabled>
>
void global_offsets(
    mpi_comm_t & comm,
    const SendType & num,
    SendType & displs,
    const IdType & recvcounts,
    const IdType & recvdispls)
{
  if (displs.size() < 1) return;
  
  auto size = displs.size() - 1;

  auto counts = make_array_ref( &displs[1], size );
  comm.all_gatherv(num, counts, recvcounts, recvdispls);
  
  displs[0] = 0;
  std::partial_sum(displs.begin(), displs.end(), displs.begin());
}
    
////////////////////////////////////////////////////////////////////////////////
/// determine offsets by exchanging counts
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
void gather_crs(
    mpi_comm_t & comm,
    const crs_t<T,U> & local,
    crs_t<T,U> & global)
{
  auto comm_size = comm.size();

  int num = local.size();
  std::vector<int> rank_counts(comm_size);
  comm.all_gather(num, rank_counts);

  std::vector<int> rank_offsets(comm_size+1);
  rank_offsets[0] = 0;
  std::partial_sum(rank_counts.begin(), rank_counts.end(), &rank_offsets[1]);

  auto local_counts = local.template counts<U>();
  
  auto & offsets = global.offsets;
  offsets.clear();
  offsets.resize(rank_offsets.back() + 1);

  auto counts = make_array_ref( &offsets[1], rank_offsets.back());
     
  comm.all_gatherv(
      local_counts, 
      counts,
      rank_counts,
      rank_offsets);

  offsets[0] = 0;
  std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());

  auto & indices = global.indices;
  indices.clear();
  indices.resize(offsets.back());

  for (int r=0; r<comm_size; ++r)
    rank_counts[r] = offsets[rank_offsets[r+1]] - offsets[rank_offsets[r]];
  std::partial_sum(rank_counts.begin(), rank_counts.end(), &rank_offsets[1]);

  comm.all_gatherv(
      local.indices, 
      indices,
      rank_counts,
      rank_offsets);
}

} // namespace

#endif
