#ifndef PUSH_QUEUE_HPP
#define PUSH_QUEUE_HPP

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
template<typename T>
class push_queue_t {
  
  mpi_comm_t comm_;
      
  const T * send_data_ = nullptr;
  const int_t * send_ids_ = nullptr;
  const int_t * send_pes_ = nullptr;
  int_t send_len_ = 0;

  std::vector<int> sendcounts_;
  std::vector<int> senddispls_;
  std::vector<int> recvcounts_;
  std::vector<int> recvdispls_;

  std::vector<T> sendbuf_;
  std::vector<T> * recvbuf_;
  
  mpi_request_t requests_;

  static bool is_async_;

  std::unique_ptr<thread_t> thread_;
  
public:

  push_queue_t(mpi_comm_t & comm);

  ~push_queue_t(){ wait(); }

  static void set_async() { 
    is_async_ = true;
  }

  void process(
      const T * data,
      const int_t * ids,
      const int_t * pes,
      int_t n,
      std::vector<T> & recvbuf);

  void wait();

  void clear();
  
  // dont call these directly unless you know what you are doing
  void _process();
  void _wait();

};

////////////////////////////////////////////////////////////////////////////////
/// Static declaraitions
////////////////////////////////////////////////////////////////////////////////
template<typename T>
bool push_queue_t<T>::is_async_ = false;

////////////////////////////////////////////////////////////////////////////////
/// Queue constructor
////////////////////////////////////////////////////////////////////////////////
template<typename T>
push_queue_t<T>::push_queue_t(mpi_comm_t & comm) :
  comm_(comm)
{  
#ifdef HAVE_THREADS
  if (is_async_) {
    thread_ = std::make_unique<thread_t>(
        [this](){ _process(); },
        [this](){ _wait(); } );
  }
#endif
}   

////////////////////////////////////////////////////////////////////////////////
/// Clear an existing queue
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void push_queue_t<T>::clear() {
  requests_.wait(); // wait for existing exchanges

  sendcounts_.clear();
  senddispls_.clear();
  recvcounts_.clear();
  recvdispls_.clear();
  sendbuf_.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// Process a batch of outgoing data
////////////////////////////////////////////////////////////////////////////////
template <typename T>
void push_queue_t<T>::process(
    const T * data,
    const int_t * ids,
    const int_t * pes,
    int_t n,
    std::vector<T> & recvbuf)
{

  send_data_ = data;
  send_ids_ = ids;
  send_pes_ = pes;
  send_len_ = n;
  recvbuf_ = &recvbuf;

#ifdef HAVE_THREADS
  if (is_async_)
    thread_->launch();
  else
    _process();
#else
  _process();
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Process a batch of outgoing data
////////////////////////////////////////////////////////////////////////////////
template <typename T>
void push_queue_t<T>::_process()
{
  clear();

  int_t comm_size = comm_.size();
  int_t comm_rank = comm_.rank();
  
  constexpr auto data_size = sizeof(T);

  // determine receiv counts
  sendcounts_.resize(comm_size, 0);
  for (int_t i=0; i<send_len_; ++i)
    sendcounts_[ send_pes_[i] ] ++;

  // exchange recvcounts
  recvcounts_.resize(comm_size);
  auto req = comm_.all_to_all(sendcounts_, recvcounts_);

  // compute displacements
  senddispls_.resize(comm_size+1);
  senddispls_[0] = 0;
  std::partial_sum(sendcounts_.begin(), sendcounts_.end(), &senddispls_[1]);

  // pack send buffer
  sendbuf_.resize(send_len_);

  std::vector<int_t> sendcounts(comm_size, 0); // temporary so we can overlap comm

  for (int_t i=0; i<send_len_; ++i) {
    auto pe = send_pes_[i];
    auto id = send_ids_[i];
    auto offset = senddispls_[pe] + sendcounts[pe];
    std::memcpy(&sendbuf_[offset], &send_data_[id], data_size);
    sendcounts[pe]++;
  }
  
  // finish recv displacements
  req.wait();

  recvdispls_.resize(comm_size+1);
  recvdispls_[0] = 0;
  std::partial_sum(recvcounts_.begin(), recvcounts_.end(), &recvdispls_[1]);

  //--- exchange
  recvbuf_->resize( recvdispls_.back() );

  for (int_t r=0; r<comm_size; ++r) {
    if (recvcounts_[r] > 0) {
      std::cout << "rank " << comm_.rank() << " receiving " << recvcounts_[r] << " from " << r 
        << " tag " << r << std::endl;
      auto ref = make_array_ref(recvbuf_->data() + recvdispls_[r], recvcounts_[r]);
      auto req = comm_.receive(ref, r, r);
      requests_.transfer(req);
    }
  }
 
  //--- send the data
  for (int_t r=0; r<comm_size; ++r) {
    if (sendcounts_[r] > 0) {
      std::cout << "rank " << comm_.rank() << " sending " << sendcounts_[r] << " to " << r 
        << " tag " << r << std::endl;
      auto ref = make_array_ref(sendbuf_.data() + senddispls_[r], sendcounts_[r]);
      auto req = comm_.send(ref, r, comm_rank);
      requests_.transfer(req);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
/// Wait for a queue to finish, and copy the results
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void push_queue_t<T>::wait() {
#ifdef HAVE_THREADS
  if (is_async_)
    thread_->wait();
  else
    _wait();
#else
  _wait();
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Wait for incoming
////////////////////////////////////////////////////////////////////////////////
template <typename T>
void push_queue_t<T>::_wait()
{
  if (requests_.empty()) return;
  
  auto is_root = comm_.is_root();

  // wait for requsts
  auto res = requests_.wait();
  if (res) {
    if (is_root)
      std::cout << "Failed waiting for communication queue, returned " << res << std::endl;
    comm_.exit(consts::failure);
  }

}


} // namespace

#endif
