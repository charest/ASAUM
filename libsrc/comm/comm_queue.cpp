#include "comm_queue.hpp"
#include "utils/mark.hpp"

#include <cstring>
#include <set>
#include <unistd.h>

namespace prl {

bool comm_queue_t::is_async_ = false;

////////////////////////////////////////////////////////////////////////////////
/// Clear an existing queue
////////////////////////////////////////////////////////////////////////////////
void comm_queue_block_t::clear() {
  requests_.wait(); // wait for existing exchanges

  data_map_.clear();
  
  sendcounts_.clear();
  senddispls_.clear();
  recvcounts_.clear();
  recvdispls_.clear();
  sendbuf_.clear();
  recvbuf_.clear();
  
  is_valid_ = false;
}

////////////////////////////////////////////////////////////////////////////////
/// process a queue
////////////////////////////////////////////////////////////////////////////////
size_t comm_queue_block_t::count(const array_ref<const int_t> & ids) const
{
  size_t sz = 0;
  for (auto & dp : data_map_.order)
    sz += dp->second->count(ids);
  return sz;
}

size_t comm_queue_block_t::pack(
    const array_ref<const int_t> & ids,
    byte_t * buf) const
{
  size_t sz = 0;
  for (auto & dp : data_map_.order)
    sz += dp->second->pack(ids, buf + sz);
  return sz;
}

size_t comm_queue_block_t::unpack(
    const array_ref<const int_t> & ids,
    const byte_t * buf) const
{
  size_t sz = 0;
  for (auto & dp : data_map_.order)
    sz += dp->second->unpack(ids, buf + sz);
  return sz;
}

void comm_queue_block_t::copy(
    const comm_queue_block_t & src_q,
    const array_ref<const int_t> & src_ids,
    const array_ref<const int_t> & dst_ids) const {
        
  auto & dst_data = data_map_.order;
  int_t num_datas = dst_data.size();
        
  const auto & src_data = src_q.data_map_.order;
    
  if (num_datas != static_cast<int_t>(src_data.size()))
    THROW_ERROR("Source and destination communication queues have "
        << "differing number of fields.");
  
  for (int_t d=0; d<num_datas; ++d) {
    auto dst_pair = dst_data[d];
    auto src_pair = src_data[d];
    dst_pair->second->copy(*src_pair->second, src_ids, dst_ids);
  } // data
}

////////////////////////////////////////////////////////////////////////////////
/// process a queue
////////////////////////////////////////////////////////////////////////////////
void comm_queue_block_t::count() {
  
  // get the shared/ghost info
  const auto & shared_ids = map_.shared_ids_;
  const auto & shared_users = map_.shared_users_;
    
  const auto & ghost_ids = map_.ghost_ids_;
  const auto & ghost_owners = map_.ghost_owners_;
    
  int_t num_shared_blocks = shared_users.size();
  int_t num_ghost_blocks = ghost_owners.size();

  // determine send counts
  sendcounts_.clear();
  sendcounts_.resize(num_shared_blocks);

  for (int_t b=0; b<num_shared_blocks; ++b)
    sendcounts_[b] = count(shared_ids.at(b));
  
  // determine receiv counts
  recvcounts_.clear();
  recvcounts_.resize(num_ghost_blocks);

  for (int_t b=0; b<num_ghost_blocks; ++b)
    recvcounts_[b] = count(ghost_ids.at(b));

  // the displacements
  int_t sendsize = sendcounts_.size();
  senddispls_.clear();
  senddispls_.resize(sendsize+1);
  senddispls_[0] = 0;
  for (int_t r=0; r<sendsize; ++r)
    senddispls_[r+1] = senddispls_[r] + sendcounts_[r];
  
  int_t recvsize = recvcounts_.size();
  recvdispls_.clear();
  recvdispls_.resize(recvsize+1);
  recvdispls_[0] = 0;
  for (int_t i=0; i<recvsize; ++i)
    recvdispls_[i+1] = recvdispls_[i] + recvcounts_[i];
  
  sendbuf_.clear();
  recvbuf_.clear();
  sendbuf_.resize(senddispls_.back());
  recvbuf_.resize(recvdispls_.back());

  is_valid_ = true;
}

////////////////////////////////////////////////////////////////////////////////
/// process a queue
////////////////////////////////////////////////////////////////////////////////
void comm_queue_block_t::process() {
  
  if (data_map_.empty()) return;

  // make sure old queue is done
  wait();

  // get the shared/ghost info
  auto shared_block_id = map_.block_id_;
  const auto & shared_ids = map_.shared_ids_;
  const auto & shared_users = map_.shared_users_;
    
  const auto & ghost_owners = map_.ghost_owners_;
  const auto & ghost_blocks = map_.ghost_blocks_;
    
  int_t num_shared_blocks = shared_users.size();
  int_t num_ghost_blocks = ghost_owners.size();

  // determine send/recv counts
  if (!is_valid_) count();
  
  //----------------------------------
  // Post recvs
  
  for (int_t b=0; b<num_ghost_blocks; ++b) {
    if (recvcounts_[b] > 0) {
      auto blk = ghost_blocks[b];
      auto rnk = ghost_owners[b];
      auto ref = make_array_ref(recvbuf_.data() + recvdispls_[b], recvcounts_[b]);
      auto req = comm_.receive(ref, rnk, blk);
      requests_.transfer(req);
    }
  }

  
  //----------------------------------
  // Pack the data
  
  auto sbuf = sendbuf_.data();

  for (int_t b=0; b<num_shared_blocks; ++b)
    pack(shared_ids.at(b), sbuf + senddispls_[b]);
  
  //----------------------------------
  // Send the data
  
  for (int_t b=0; b<num_shared_blocks; ++b) {
    if (sendcounts_[b] > 0) {
      auto ref = make_array_ref(sendbuf_.data() + senddispls_[b], sendcounts_[b]);
      auto r = shared_users[b];
      auto req = comm_.send(ref, r, shared_block_id);
      requests_.transfer(req);
    }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
/// Wait for a queue to finish, and copy the results
////////////////////////////////////////////////////////////////////////////////
void comm_queue_block_t::wait() {

  if (data_map_.empty()) return;
  if (requests_.empty()) return;

  // wait for requsts
  auto res = requests_.wait();
  if (res) {
    if (comm_.is_root()) {
      std::cout << "Failed waiting for communication queue, returned " << res << std::endl;
    }
    comm_.exit(consts::failure);
  }
  
  // get ghost info
  const auto & ghost_ids = map_.ghost_ids_;
  int_t num_ghost_blocks = ghost_ids.size();
  
  // Unpack the data
  auto buf = recvbuf_.data();

  for (int_t b=0; b<num_ghost_blocks; ++b)
    unpack(ghost_ids.at(b), buf + recvdispls_[b]);
  
}

////////////////////////////////////////////////////////////////////////////////
/// Queue constructor
////////////////////////////////////////////////////////////////////////////////
comm_queue_t::comm_queue_t(mpi_comm_t & comm, comm_map_t & comm_maps) :
  parent_comm_(comm)
{  

#ifdef HAVE_THREADS
  if (is_async_) {
    thread_ = std::make_unique<thread_t>(
        [this](){ _process(); },
        [this](){ _wait();  } );
  }
#endif

  queue_index_.resize(comm_maps.size());
  int_t i=0;
  for (auto & cm : comm_maps) {
    auto it = queues_.emplace( queues_.end(), parent_comm_, cm);
    queue_map_[cm.block_id_] = it; 
    queue_index_[i] = it;
    i++;
  }

}   

////////////////////////////////////////////////////////////////////////////////
/// Clear an existing queue
////////////////////////////////////////////////////////////////////////////////
void comm_queue_t::clear() {
  request_.wait(); // wait for existing exchanges

  sendcounts_.clear();
  senddispls_.clear();
  recvcounts_.clear();
  recvdispls_.clear();
  sendbuf_.clear();
  recvbuf_.clear();

  comm_.reset();

  is_valid_ = false;
  
  for (auto & q : queues_) q.clear();
}

////////////////////////////////////////////////////////////////////////////////
/// process a queue
////////////////////////////////////////////////////////////////////////////////
void comm_queue_t::process() {
  MARK_FUNCTION();
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
/// Compute send/recv counts
////////////////////////////////////////////////////////////////////////////////
void comm_queue_t::count() {
  
  int_t world_rank = parent_comm_.rank();

  // determine the ranks with data
  comm_ = std::make_unique<sub_comm_t>(parent_comm_, queues_.size());
  int_t comm_size = comm_->comm().size();

  // determine receiv counts
  recvcounts_.clear();
  recvcounts_.resize(comm_size, 0);

  for (auto & q : queues_) {
    const auto & map  = q.map_;
    const auto & ghost_ids = map.ghost_ids_;
    const auto & ghost_owners = map.ghost_owners_;
    const auto & ghost_blocks = map.ghost_blocks_;
    int_t num_blks = ghost_owners.size();
    for (int_t b=0; b<num_blks; ++b) {
      auto r = ghost_owners[b];
      if (r != world_rank) {
        auto sz = q.count(ghost_ids.at(b)); 
        auto r2 = comm_->map_rank(r);
        auto bid = ghost_blocks[b];
        if (r2 < 0) THROW_ERROR("Rank maps to one with no data!");
        recvcounts_[r2] += sz + 2*sizeof(decltype(bid));
      }
    } // block
  } // queue


  // determine send counts
  sendcounts_.clear();
  sendcounts_.resize(comm_size, 0);

  for (auto & q : queues_) {
    const auto & map  = q.map_;
    const auto & shared_ids = map.shared_ids_;
    const auto & shared_users = map.shared_users_;
    const auto & shared_blocks = map.shared_blocks_;
    int_t num_shared_blocks = shared_ids.size();
    for (int_t b=0; b<num_shared_blocks; ++b) {
      auto r = shared_users[b];
      if (r != world_rank) {
        auto sz = q.count(shared_ids.at(b)); 
        auto i = comm_->map_rank(r);
        auto bid = shared_blocks[b];
        if (i < 0) THROW_ERROR("Rank maps to one with no data!");
        sendcounts_[i] += sz + 2*sizeof(decltype(bid));
      }
    } // block 
  } // queues

  // the displacements
  senddispls_.clear();
  recvdispls_.clear();
  senddispls_.resize(comm_size+1);
  recvdispls_.resize(comm_size+1);
  senddispls_[0] = 0;
  recvdispls_[0] = 0;
  for (int r=0; r<comm_size; ++r) {
    senddispls_[r+1] = senddispls_[r] + sendcounts_[r];
    recvdispls_[r+1] = recvdispls_[r] + recvcounts_[r];
  }
  
  sendbuf_.clear();
  recvbuf_.clear();
  sendbuf_.resize(senddispls_[comm_size]);
  recvbuf_.resize(recvdispls_[comm_size]);

  is_valid_ = true;
}

////////////////////////////////////////////////////////////////////////////////
/// process a queue
////////////////////////////////////////////////////////////////////////////////
void comm_queue_t::_process() {
  
  // make sure old queue is done
  _wait();

  int_t world_rank = parent_comm_.rank();
  
  //----------------------------------
  // determine send/recv counts
  if (!is_valid_) count();


  // empty ranks can leave
  if (queues_.empty()) return;
  
  //----------------------------------
  // Pack data
  
  auto buf = sendbuf_.data();
  std::fill(sendcounts_.begin(), sendcounts_.end(), 0);

  for (auto & q : queues_) {
    const auto & map  = q.map_;
    const auto & shared_ids = map.shared_ids_;
    const auto & shared_users = map.shared_users_;
    const auto & shared_blocks = map.shared_blocks_;
    auto shared_block_id = map.block_id_;
    int_t num_shared_blocks = shared_ids.size();
  
    for (int_t b=0; b<num_shared_blocks; ++b) {  
      auto r = shared_users[b];
      if (r != world_rank) {
        auto r2 = comm_->map_rank(r);
        if (r2 < 0) THROW_ERROR("Rank maps to one with no data!");
        
        auto dst1 = buf + senddispls_[r2] + sendcounts_[r2];
        auto sz = sizeof(decltype(shared_block_id));
        std::memcpy(dst1, &shared_block_id, sz);
        sendcounts_[r2] += sz;
            
        auto dst2 = buf + senddispls_[r2] + sendcounts_[r2];
        auto dst_block_id = shared_blocks[b];
        sz = sizeof(decltype(dst_block_id));
        std::memcpy(dst2, &dst_block_id, sz);
        sendcounts_[r2] += sz;

        sendcounts_[r2] += q.pack(
            shared_ids.at(b),
            buf + senddispls_[r2] + sendcounts_[r2]
        );

      } // rank
    } // blocks

  } // queue

  //----------------------------------
  // Send the data
 
  auto res = comm_->comm().all_to_allv(
      sendbuf_,
      sendcounts_,
      senddispls_,
      recvbuf_,
      recvcounts_,
      recvdispls_);
  request_ = std::move(res);
  
  //----------------------------------
  // Copy ranks own data

  for (auto & dst_q : queues_) {
    const auto & dst_map  = dst_q.map_;
    auto dst_block_id = dst_map.block_id_;
    const auto & ghost_owners = dst_map.ghost_owners_;
    const auto & ghost_blocks = dst_map.ghost_blocks_;
    const auto & ghost_ids = dst_map.ghost_ids_;
    int_t num_blks = ghost_owners.size();


    for (int_t bdst=0; bdst<num_blks; ++bdst) {
      auto r = ghost_owners[bdst];
      if (r == world_rank) {
  
        auto src_block_id = ghost_blocks[bdst];
        const auto & src_q = *queue_map_.at(src_block_id);

        const auto & src_map = src_q.map_;
        const auto & shared_ids = src_map.shared_ids_;
        auto bsrc = src_map.shared_block_map_.at(dst_block_id);

        auto sids = shared_ids.at(bsrc);
        auto num_shared = sids.size();

        auto gids = ghost_ids.at(bdst);
        auto num_ghost = gids.size();

        if (num_ghost != num_shared) {
          THROW_ERROR("Ghost and shared sizes do not match.  "
            << "block " << dst_block_id << " has " << num_ghost
            << " ghost but block " << src_block_id << " has "
            << num_shared);
        }

        dst_q.copy(src_q, sids, gids);

      } // rank

    } // block
    
  } // queues
}

////////////////////////////////////////////////////////////////////////////////
/// Wait for a queue to finish, and copy the results
////////////////////////////////////////////////////////////////////////////////
void comm_queue_t::wait() {
  MARK_FUNCTION();
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
/// Wait for a queue to finish, and copy the results
////////////////////////////////////////////////////////////////////////////////
void comm_queue_t::_wait() {
  if (request_.empty()) return;

  // wait for requsts
  auto res = request_.wait();
  if (res) {
    if (comm_->comm().is_root()) {
      std::cout << "Failed waiting for communication queue, returned " << res << std::endl;
    }
    comm_->comm().exit(consts::failure);
  }

  int_t comm_size = comm_->comm().size();
  auto buf = recvbuf_.data();
    
  for (int_t r=0; r<comm_size; ++r) { 
    auto cnt = recvdispls_[r];
    for (; cnt<recvdispls_[r+1];) {

      // get source and destination block id
      int_t src_block_id;
      std::memcpy(&src_block_id, &buf[cnt], sizeof(int_t));
      cnt += sizeof(int_t);
      
      int_t dst_block_id;
      std::memcpy(&dst_block_id, &buf[cnt], sizeof(int_t));
      cnt += sizeof(int_t);
    
      // find the correct queue
      const auto & q = *queue_map_.at(dst_block_id);
      const auto & map  = q.map_;
      const auto & ghost_owners = map.ghost_owners_;
      const auto & ghost_ids = map.ghost_ids_;
      int_t num_blks = ghost_owners.size();

      // and its destination block
      auto b = map.ghost_block_map_.at(src_block_id);
      if (b >= num_blks) {
        THROW_ERROR("The desitnation block is out of bounds.  "
            << "Block " << dst_block_id << " only has "
            << num_blks << " ghost blocks, you are asking for the "
            << b << "th block.");
      }
  
      // Unpack the data
      auto rnk = ghost_owners[b];
      if (rnk != r) {
        THROW_ERROR("The ghost block came from rank "
            << r << " but is supposed to have come from "
            << rnk << ".");
      }
        
      auto r2 = comm_->map_rank(r);
      if (r2 < 0) THROW_ERROR("Rank maps to one with no data!");

      cnt += q.unpack(ghost_ids.at(b), buf + cnt);

    } // displs

    if (cnt != recvdispls_[r+1]) {
      THROW_ERROR("Receive counts do not match the expected values.");
    }

  } // comm
  
  
}

} // namespace
