#ifndef DIRECTED_QUEUE_HPP
#define DIRECTED_QUEUE_HPP

#include "config.hpp"
#include "comm_map.hpp"
#include "mpi_comm.hpp"
#include "sub_comm.hpp"

#include "utils/ordered_map.hpp"
#include "utils/thread.hpp"

#include <cstring>
#include <map>
#include <memory>
#include <vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// The communication queue
////////////////////////////////////////////////////////////////////////////////
struct directed_queue_msg_info_t {
  int_t tag;
  int_t rank;
  mpi_request_t request;
  std::vector<byte_t> buffer;
};

////////////////////////////////////////////////////////////////////////////////
/// The communication queue
////////////////////////////////////////////////////////////////////////////////
class directed_queue_t {
  
  mpi_comm_t comm_;

  std::vector<directed_queue_msg_info_t> received_messages_;
  std::vector<directed_queue_msg_info_t> sent_messages_;
      
  //const T * send_data_ = nullptr;
  const int_t * send_ids_ = nullptr;
  const int_t * send_pes_ = nullptr;
  int_t send_len_ = 0;

  std::vector<int> sendcounts_;
  std::vector<int> senddispls_;
  std::vector<int> recvcounts_;
  std::vector<int> recvdispls_;

  //std::vector<T> sendbuf_;
  //std::vector<T> * recvbuf_;
  
  mpi_request_t requests_;

  static bool is_async_;

  thread_t thread_;
  
public:

  directed_queue_t(mpi_comm_t & comm);

  ~directed_queue_t(){ wait(); }

  static void set_async() { 
    is_async_ = true;
  }

#if 0
  void process(
      const T * data,
      const int_t * ids,
      const int_t * pes,
      int_t n,
      std::vector<T> & recvbuf);
#endif

  template<typename T>
  void send(const T & data, int_t rank, int_t tag)
  {
    directed_queue_msg_info_t info;
    info.tag = tag;
    info.rank = rank;
    auto data_size = data.size() * sizeof(T);
    info.buffer.resize(data_size);
    std::memcpy(info.buffer.data(), data.data(), data_size);
    info.request = comm_.send(info.buffer, rank, tag);
  }

  const std::vector<directed_queue_msg_info_t> & messages() const
  { return received_messages_; }

  void wait();

  void clear();
  
  // dont call these directly unless you know what you are doing
  void _process();
  void _wait();

};

} // namespace

#endif
