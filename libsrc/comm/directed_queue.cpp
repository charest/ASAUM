#include "directed_queue.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Static declaraitions
////////////////////////////////////////////////////////////////////////////////
bool directed_queue_t::is_async_ = false;

////////////////////////////////////////////////////////////////////////////////
/// Queue constructor
////////////////////////////////////////////////////////////////////////////////
directed_queue_t::directed_queue_t(mpi_comm_t & comm) :
  comm_(comm),
  thread_( [this](){ _process(); }, [this](){ _wait(); } )
{  
}   

////////////////////////////////////////////////////////////////////////////////
/// Clear an existing queue
////////////////////////////////////////////////////////////////////////////////
void directed_queue_t::clear() {
  requests_.wait(); // wait for existing exchanges

  sendcounts_.clear();
  senddispls_.clear();
  recvcounts_.clear();
  recvdispls_.clear();
  //sendbuf_.clear();
}

#if 0
////////////////////////////////////////////////////////////////////////////////
/// Process a batch of outgoing data
////////////////////////////////////////////////////////////////////////////////
void directed_queue_t::process(
    const v * data,
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
    thread_.launch();
  else
    _process();
#else
  _process();
#endif
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// Process a batch of outgoing data
////////////////////////////////////////////////////////////////////////////////
void directed_queue_t::_process()
{
#if 0
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
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Wait for a queue to finish, and copy the results
////////////////////////////////////////////////////////////////////////////////
void directed_queue_t::wait() {
#ifdef HAVE_THREADS
  if (is_async_)
    thread_.wait();
  else
    _wait();
#else
  _wait();
#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Wait for incoming
////////////////////////////////////////////////////////////////////////////////
void directed_queue_t::_wait()
{
#if 0
  if (requests_.empty()) return;
  
  auto is_root = comm_.is_root();

  // wait for requsts
  auto res = requests_.wait();
  if (res) {
    if (is_root)
      std::cout << "Failed waiting for communication queue, returned " << res << std::endl;
    comm_.exit(consts::failure);
  }
#endif
}

} // namespace
