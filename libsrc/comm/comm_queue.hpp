#ifndef COMM_QUEUE_HPP
#define COMM_QUEUE_HPP

#include "config.hpp"
#include "comm_info.hpp"
#include "comm_map.hpp"
#include "mpi_comm.hpp"
#include "sub_comm.hpp"

#include "utils/ordered_map.hpp"
#include "utils/thread.hpp"

#include <map>
#include <memory>
#include <vector>

namespace prl {

class comm_queue_t;
  
////////////////////////////////////////////////////////////////////////////////
/// The block level communication queue
////////////////////////////////////////////////////////////////////////////////
class comm_queue_block_t {

  using map_type = ordered_map<void *, std::unique_ptr<comm_info_t>>;
  map_type data_map_;

  std::vector<int> sendcounts_;
  std::vector<int> senddispls_;
  std::vector<int> recvcounts_;
  std::vector<int> recvdispls_;

  std::vector<byte_t> sendbuf_;
  std::vector<byte_t> recvbuf_;

  mpi_comm_t & comm_;
  const comm_map_block_t & map_;
 
  mpi_request_t requests_;
  
  bool is_valid_ = false;

  void reset() { is_valid_ = false; }

  void count();
  size_t count(const array_ref<const int_t> & ids) const;
  size_t pack(const array_ref<const int_t> & ids, byte_t * buf) const;
  size_t unpack(const array_ref<const int_t> & ids, const byte_t * buf) const;

  void copy(
    const comm_queue_block_t & src_q,
    const array_ref<const int_t> & src_ids,
    const array_ref<const int_t> & dst_ids) const;

public:

  comm_queue_block_t( mpi_comm_t & comm, const comm_map_block_t & map) :
    comm_(comm), map_(map)
  {}

  ~comm_queue_block_t() { wait(); }
  
  template<typename T>
  void add( T & data )
  {
    if (data.empty()) return;
    reset();
    auto ptr = static_cast<void*>(data.data());
    if (data.size() % map_.index_size_ != 0)
      THROW_ERROR("Data added to comomunication queue does not divide evenly "
          << " by the index size" );
    auto data_size = data.size() / map_.index_size_ * sizeof(typename T::value_type);
    data_map_.emplace( ptr, std::make_unique<dense_info_t>(data_size, ptr) );
  }
  
  template<typename T, typename U>
  void add(T & data, U & offsets)
  {
    if (offsets.empty()) {
      add(data);
      return;
    }

    reset();
    auto ptr = static_cast<void*>(data.data());
    auto data_size = offsets.size() - 1;
    if (data_size % map_.index_size_ != 0)
      THROW_ERROR("Data added to comomunication queue does not divide evenly "
          << " by the index size" );
    data_size = data_size / map_.index_size_ * sizeof(typename T::value_type);
    data_map_.emplace(
        ptr,
        std::make_unique<crs_info_t>(data_size, ptr, offsets.data())
    );
  }
  
  template<typename T, typename U, typename V>
  void add(T & data, U & offsets, V & counts)
  {
    if (offsets.empty()) {
      add(data);
      return;
    }

    if (counts.empty()) {
      add(data, offsets);
      return;
    }

    reset();
    auto ptr = static_cast<void*>(data.data());
    if (counts.size() % map_.index_size_ != 0)
      THROW_ERROR("Data added to comomunication queue does not divide evenly "
          << " by the index size" );
    auto data_size = counts.size() / map_.index_size_ * sizeof(typename T::value_type);
    data_map_.emplace(
        ptr,
        std::make_unique<crs_info_t>(data_size, ptr, offsets.data(), counts.data())
    );
  }


  void process();

  void wait();

  void clear();

  friend comm_queue_t;

};

////////////////////////////////////////////////////////////////////////////////
/// The communication queue
////////////////////////////////////////////////////////////////////////////////
class comm_queue_t {
  
  using list_type = std::list<comm_queue_block_t>;

  mpi_comm_t parent_comm_;
  std::unique_ptr<sub_comm_t> comm_;

  list_type queues_;
  std::vector<list_type::iterator> queue_index_;
  std::map<int_t, list_type::iterator> queue_map_;
  
  std::vector<int> sendcounts_;
  std::vector<int> senddispls_;
  std::vector<int> recvcounts_;
  std::vector<int> recvdispls_;

  std::vector<byte_t> sendbuf_;
  std::vector<byte_t> recvbuf_;
  
  mpi_request_t request_;
  
  bool is_valid_ = false;

  static bool is_async_;

  std::unique_ptr<thread_t> thread_;
  
  void reset() {
    is_valid_ = false;
    for (auto & q : queues_) q.reset();
  }
  
  void count();

public:

  comm_queue_t(mpi_comm_t & comm, comm_map_t & comm_maps);
  
  ~comm_queue_t() { wait(); }

  comm_queue_block_t & operator[](int_t i) 
  { return *queue_index_[i]; }

  static void set_async() { is_async_ = true; }

  auto begin() { return queues_.begin(); }
  auto end() { return queues_.end(); }

  void process();
  void wait();

  void clear();
  
  // dont call these directly unless you know what you are doing
  void _process();
  void _wait();

};


} // namespace

#endif
