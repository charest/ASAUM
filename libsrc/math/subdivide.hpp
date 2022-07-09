#ifndef SUBDIVIDE_HPP
#define SUBDIVIDE_HPP

#include "comm/mpi_comm.hpp"
#include "utils/string_utils.hpp"

#include<algorithm>
#include<numeric>
#include<map>
#include<vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// \brief A simple utility for subdividing an index space into several parts 
////////////////////////////////////////////////////////////////////////////////
template< typename T, typename U >
void subdivide( T nelem, U npart, std::vector<T> & dist ) {

  T quot = nelem / npart;
  T rem = nelem % npart;

  dist.reserve( npart+1 );
  dist.push_back(0);

  // Set the distributions for each rank. This happens on all ranks.
  // Each rank gets the average number of indices, with higher ranks
  // getting an additional index for non-zero remainders.
  for(U r(0); r < npart; ++r) {
    T indices = quot + ((static_cast<T>(r) < rem) ? 1 : 0);
    dist.push_back(dist[r] + indices);
  } // for
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Simple utility for determining which rank owns an id
////////////////////////////////////////////////////////////////////////////////
template<typename T>
T owner( const std::vector<T> & distribution, T i )
{
  auto it = std::upper_bound( distribution.begin(), distribution.end(), i );
  return std::distance( distribution.begin(), it ) - 1;
}
  
template<typename T>
void invert_owner(const std::vector<T> & dist, std::vector<T> & out)
{
  out.clear();
  T n = dist.size();
  if (n==0) return;

  auto l = dist.back();
  out.resize(l, -1);
  
  for (T i=0, r=0; i<l && r<n-1; ++i) {
    while (i>=dist[r+1]) ++r;
    out[i] = r;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// is an id owned by this rank
////////////////////////////////////////////////////////////////////////////////
bool owned_by( const std::vector<size_t> & dist, size_t id, size_t owner);


////////////////////////////////////////////////////////////////////////////////
/// \brief Simple utility for determining which rank owns an id
/// \remark sequential distribution
////////////////////////////////////////////////////////////////////////////////
template<typename T>
std::vector<T> distribute_sequentially(
    T nelem,
    T npart,
    T ipart)
{
  std::vector<T> results;
  
  if (npart >= nelem) {
    if (ipart < nelem)
      results.emplace_back(ipart);
  }
  else if (ipart<npart) {
    T quot = nelem / npart;
    T rem = nelem % npart;
    T num = quot + ((ipart<rem) ? 1 : 0);
    T istart = ipart*quot + std::min(ipart, rem);
    results.resize(num);
    std::iota(results.begin(), results.end(), istart);
  }
  return results;
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Simple utility for determining which rank owns an id
/// \remark roundrobin distribution
////////////////////////////////////////////////////////////////////////////////
template<typename T>
std::vector<T> distribute_cyclic(
    T nelem,
    T npart,
    T ipart)
{
  std::vector<T> results;
  
  if (npart >= nelem) {
    auto n = npart / nelem;
    auto r = ipart / n;
    auto q = ipart % n;
    if (q == 0 && r < nelem)
      results.emplace_back(r);
  }
  else if (ipart<npart) {
    T quot = nelem / npart;
    T rem = nelem % npart;
    T num = quot + ((ipart<rem) ? 1 : 0);
    results.resize(num);
    results[0] = ipart;
    for (T i=0; i<num-1; ++i)
      results[i+1] = results[i] + quot;
  }
  return results;
}

////////////////////////////////////////////////////////////////////////////////
/// \brief Simple utility for determining which rank owns an id
/// \remark based on a key/hostname distribution
////////////////////////////////////////////////////////////////////////////////
template<typename T>
std::vector<T> distribute_key(
    T nelem,
    mpi_comm_t & comm,
    const std::string & key)
{
  std::vector<T> results;
  T ipart = comm.rank();
  T npart = comm.size();
 
  auto hsh = string_hash<size_t>(key, key.length());
    
  std::vector<size_t> hashes(npart);
  auto req = comm.all_gather(hsh, hashes);
  req.wait();

  std::map<size_t, T> hostmap;
  for (T i=0; i<npart; ++i)
    hostmap.emplace( hashes[i], hostmap.size() );
  
  T num_hosts = hostmap.size();

  std::vector<T> local_rank_order(npart);
  std::vector<T> host_counts(num_hosts, 0);
  for (T i=0; i<npart; ++i) {
    auto j = hostmap.at( hashes[i] );
    local_rank_order[i] = host_counts[j];
    host_counts[j]++;
  }

  auto n = nelem / num_hosts;
  auto q = nelem % num_hosts;
  
  auto host_id = hostmap.at(hsh);
  auto host_start = n * host_id;
  host_start += std::min<size_t>(q, host_id);
  
  T num_local = host_id < q ? n+1 : n;
  
  if (npart >= nelem) {
    auto local_id = local_rank_order[ipart];
    if (local_id<num_local)
      results.emplace_back(host_start + local_id);
  }
  else if (ipart<npart) {
    std::vector<T> local_div;
    subdivide(num_local, host_counts[host_id], local_div);
    auto local_id = local_rank_order[ipart];
    auto num = local_div[local_id+1] - local_div[local_id];
    auto start = host_start + local_div[local_id];
    results.resize(num);
    std::iota(results.begin(), results.end(), start);
  }
  return results;
}


} // namespace

#endif
