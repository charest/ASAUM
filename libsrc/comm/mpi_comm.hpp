#ifndef MPI_COMMUNICATOR_HPP
#define MPI_COMMUNICATOR_HPP

#include "comm_mutex.hpp"
#include "mpi_types.hpp"
#include "utils/type_traits.hpp"

#include <functional>
#include <vector>

namespace prl {
  
////////////////////////////////////////////////////////////////////////////////
/// Reduction operations
////////////////////////////////////////////////////////////////////////////////
enum class redop_t {
  max,
  min,
  sum,
  prod,
  lor,
  land,
  bor,
  band
};

/// the mpi datatype
using mpi_data_t = MPI_Datatype*;

////////////////////////////////////////////////////////////////////////////////
/// MPI request type
////////////////////////////////////////////////////////////////////////////////
class mpi_request_t {
  std::vector<MPI_Request> requests_;

  MPI_Request * create();

public:

  mpi_request_t() = default;
  mpi_request_t(const mpi_request_t&) = delete;
  mpi_request_t(mpi_request_t&&) = default;
  mpi_request_t & operator=(mpi_request_t&&) = default;

  void transfer(mpi_request_t & req);

  int wait();

  bool empty() { return requests_.empty(); }

  bool finished();

  ~mpi_request_t() { wait(); }
  
  friend class mpi_comm_t;
};
      

////////////////////////////////////////////////////////////////////////////////
/// MPI function
////////////////////////////////////////////////////////////////////////////////
class mpi_function_t {

  typedef void (*fun_type)(void *, void *, int *, mpi_data_t);
  std::shared_ptr<MPI_Op> op_;

  MPI_Op * create();

public:

  mpi_function_t() = default;

  MPI_Op get() const { 
    if (op_) return *op_;
    else return MPI_OP_NULL;
  }
  
  friend class mpi_comm_t;
};

////////////////////////////////////////////////////////////////////////////////
/// MPI communicator
////////////////////////////////////////////////////////////////////////////////
class mpi_comm_t {

private:

  MPI_Comm comm_ = MPI_COMM_NULL;
  MPI_Group group_ = MPI_GROUP_NULL;
  int rank_ = -1;
  int size_ = -1;
    
  char hostname_[MPI_MAX_PROCESSOR_NAME];

  bool top_ = false;
  
  static MPI_Op redop_id(redop_t op) {
    switch (op) {
      case (redop_t::max):  return  MPI_MAX;
      case (redop_t::min):  return  MPI_MIN;
      case (redop_t::sum):  return  MPI_SUM;
      case (redop_t::prod): return  MPI_PROD;
      case (redop_t::lor):  return  MPI_LOR;
      case (redop_t::land): return  MPI_LAND;
      case (redop_t::bor):  return  MPI_BOR;
      case (redop_t::band): return  MPI_BAND;
    };
    return MPI_OP_NULL;
  }

  void shutdown();

public:

  mpi_comm_t(int argc, char * argv[], bool is_threaded = false);
  
  mpi_comm_t(const mpi_comm_t & comm);

  mpi_comm_t(const mpi_comm_t & comm, bool flag);

  mpi_comm_t(
    const mpi_comm_t & comm,
    const std::vector<char> & rank_flags);

  //mpi_comm_t(const mpi_comm_t &) = delete;

  virtual ~mpi_comm_t();

  int rank() const { return rank_; }
  int size() const { return size_; }
  MPI_Comm comm() const { return comm_; }

  const char * hostname() const { return hostname_; }

  bool is_root() const { return rank_ == 0; }
  
  void exit(int i);

  void barrier();

  bool probe(int source, int tag);
  bool probe();
  bool probe_any(int & source, int & tag);
  bool probe_tag(int & source, int tag);

  mpi_function_t create_op(
      mpi_function_t::fun_type fun,
      bool commutes = true );

  template<typename T>
  mpi_request_t all_reduce(const T & in, T & out, redop_t op)
  {
    auto opid = redop_id(op);
    auto type_id = mpi_type_t<T>::value();
    
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Iallreduce(&in, &out, 1, type_id, opid, comm_, req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }
  
  template<
    typename InType,
    typename OutType,
    typename ValueType = typename InType::value_type,
    bool Enabled =
      is_container_v<InType> &&
      is_container_v<OutType> &&
      std::is_same<ValueType, typename OutType::value_type>::value,
    typename = typename std::enable_if_t<Enabled>
  >
  mpi_request_t all_reduce(
      const InType & in,
      OutType & out,
      redop_t op)
  {
    auto opid = redop_id(op);
    auto type_id = mpi_type_t<ValueType>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Iallreduce(in.data(), out.data(), in.size(), type_id, opid, comm_, req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }
  
  template<typename T>
  mpi_request_t all_reduce(const T & in, T & out, const mpi_function_t & op)
  {
    auto n = sizeof(T);
    auto opid = op.get();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Iallreduce(&in, &out, n, MPI_BYTE, opid, comm_, req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }
    
  template<
    typename InType,
    typename OutType,
    bool Enabled =
      is_container_v<InType> &&
      is_container_v<OutType>,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t all_to_all(
      const InType & in,
      OutType & out,
      int cnt = 1)
  {
    using send_type = typename InType::value_type;
    using recv_type = typename OutType::value_type;
    auto in_type = mpi_type_t<send_type>::value();
    auto out_type = mpi_type_t<recv_type>::value();
    auto in_cnt = cnt;
    auto out_cnt = cnt * sizeof(send_type) / sizeof(recv_type);
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Ialltoall( in.data(), in_cnt, in_type, out.data(),
        out_cnt, out_type, comm_, req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }

  template<
    typename SendType,
    typename RecvType,
    typename SendIdType,
    typename RecvIdType,
    bool Enabled =
      is_container_v<SendType> &&
      is_container_v<RecvType> &&
      is_container_v<SendIdType> &&
      is_container_v<RecvIdType> &&
      std::is_same<int, typename SendIdType::value_type>::value &&
      std::is_same<int, typename RecvIdType::value_type>::value,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t all_to_allv(
      const SendType & sendbuf,
      const SendIdType & sendcounts,
      const SendIdType & senddispls,
      RecvType & recvbuf,
      const RecvIdType & recvcounts,
      const RecvIdType & recvdispls)
  {
    auto sendtype = mpi_type_t<typename SendType::value_type>::value();
    auto recvtype = mpi_type_t<typename RecvType::value_type>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Ialltoallv(
        sendbuf.data(),
        sendcounts.data(),
        senddispls.data(),
        sendtype,
        recvbuf.data(),
        recvcounts.data(),
        recvdispls.data(),
        recvtype,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }

  template<
    typename T,
    bool Enabled = !is_container_v<T>,
    typename std::enable_if_t<Enabled>* = nullptr>
  mpi_request_t send(
      const T & sendbuf,
      int dest,
      int tag = 0)
  {
    auto type = mpi_type_t<T>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Isend(
        &sendbuf,
        1,
        type,
        dest,
        tag,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }

  template<
    typename T,
    bool Enabled = is_container_v<T>,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t send(
      const T & sendbuf,
      int dest,
      int tag = 0)
  {
    auto type = mpi_type_t<typename T::value_type>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Isend(
        sendbuf.data(),
        sendbuf.size(),
        type,
        dest,
        tag,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }

  template<
    typename T,
    bool Enabled = !is_container_v<T>,
    typename std::enable_if_t<Enabled>* = nullptr>
  mpi_request_t receive(
      T & recvbuf,
      int src,
      int tag = MPI_ANY_TAG)
  {
    auto type = mpi_type_t<T>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Irecv(
        &recvbuf,
        1,
        type,
        src,
        tag,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }

  template<
    typename T,
    bool Enabled = is_container_v<T>,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t receive(
      T & recvbuf,
      int src,
      int tag = MPI_ANY_TAG)
  {
    auto type = mpi_type_t<typename T::value_type>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Irecv(
        recvbuf.data(),
        recvbuf.size(),
        type,
        src,
        tag,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }
  
  template<
    typename SendType,
    typename RecvType,
    bool Enabled =
      !is_container_v<SendType> &&
      is_container_v<RecvType>,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t gather(
      const SendType & sendbuf,
      RecvType & recvbuf,
      int dst)
  {
    auto sendtype = mpi_type_t<SendType>::value();
    auto recvtype = mpi_type_t<typename RecvType::value_type>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Igather(
        &sendbuf,
        1,
        sendtype,
        recvbuf.data(),
        1,
        recvtype,
        dst,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }
  
  template<
    typename T,
    bool Enabled = !is_container_v<T>,
    typename std::enable_if_t<Enabled>* = nullptr>
  mpi_request_t gather(
      const T & sendbuf,
      int dst)
  {
    auto sendtype = mpi_type_t<T>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Igather(
        &sendbuf,
        1,
        sendtype,
        nullptr,
        0,
        sendtype,
        dst,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }


  
  template<
    typename SendType,
    typename RecvType,
    typename IdType,
    bool Enabled =
      is_container_v<SendType> &&
      is_container_v<RecvType> &&
      is_container_v<IdType> &&
      std::is_same<int, typename IdType::value_type>::value,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t gatherv(
      const SendType & sendbuf,
      RecvType & recvbuf,
      const IdType & recvcounts,
      const IdType & recvdispls,
      int dst)
  {
    auto sendtype = mpi_type_t<typename SendType::value_type>::value();
    auto recvtype = mpi_type_t<typename RecvType::value_type>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Igatherv(
        sendbuf.data(),
        sendbuf.size(),
        sendtype,
        recvbuf.data(),
        recvcounts.data(),
        recvdispls.data(),
        recvtype,
        dst,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }

  template<
    typename T,
    bool Enabled = is_container_v<T>,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t gatherv(
      const T & sendbuf,
      int dst)
  {
    auto sendtype = mpi_type_t<typename T::value_type>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Igatherv(
        sendbuf.data(),
        sendbuf.size(),
        sendtype,
        nullptr,
        nullptr,
        nullptr,
        sendtype,
        dst,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }


  template<
    typename InType,
    typename OutType,
    bool Enabled = is_container_v<OutType>,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t all_gather(const InType & in, OutType & out)
  {
    using recv_type = typename OutType::value_type;
    auto in_type = mpi_type_t<InType>::value();
    auto out_type = mpi_type_t<recv_type>::value();
    auto out_cnt = sizeof(InType) / sizeof(recv_type);
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Iallgather( &in, 1, in_type, out.data(),
        out_cnt, out_type, comm_, req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }
  
  template<
    typename SendType,
    typename RecvType,
    typename IdType,
    bool Enabled =
      !is_container_v<SendType> &&
      is_container_v<RecvType> &&
      is_container_v<IdType> &&
      std::is_same<int, typename IdType::value_type>::value,
    typename std::enable_if_t<Enabled>* = nullptr>
  mpi_request_t all_gatherv(
      const SendType & sendbuf,
      RecvType & recvbuf,
      const IdType & recvcounts,
      const IdType & recvdispls)
  {
    auto sendtype = mpi_type_t<SendType>::value();
    auto recvtype = mpi_type_t<typename RecvType::value_type>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Iallgatherv(
        &sendbuf,
        1,
        sendtype,
        recvbuf.data(),
        recvcounts.data(),
        recvdispls.data(),
        recvtype,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }

  
  template<
    typename SendType,
    typename RecvType,
    typename IdType,
    bool Enabled =
      is_container_v<SendType> &&
      is_container_v<RecvType> &&
      is_container_v<IdType> &&
      std::is_same<int, typename IdType::value_type>::value,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t all_gatherv(
      const SendType & sendbuf,
      RecvType & recvbuf,
      const IdType & recvcounts,
      const IdType & recvdispls)
  {
    auto sendtype = mpi_type_t<typename SendType::value_type>::value();
    auto recvtype = mpi_type_t<typename RecvType::value_type>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Iallgatherv(
        sendbuf.data(),
        sendbuf.size(),
        sendtype,
        recvbuf.data(),
        recvcounts.data(),
        recvdispls.data(),
        recvtype,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }

  
  template<
    typename T,
    bool Enabled = !is_container_v<T>,
    typename std::enable_if_t<Enabled>* = nullptr>
  mpi_request_t broadcast(T & buffer, int src)
  {
    auto type = mpi_type_t<T>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Ibcast(&buffer, 1, type, src, comm_, req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    
    return req;
  }
  
  template<
    typename T,
    bool Enabled = is_container_v<T>,
    typename = typename std::enable_if_t<Enabled>>
  mpi_request_t broadcast(T & buffer, int src)
  {
    using value_type = typename T::value_type;
    auto type = mpi_type_t<value_type>::value();
    mpi_request_t req;
    comm_mutex_t::instance().lock();
    auto ret = MPI_Ibcast(
        const_cast<value_type*>(buffer.data()), // for std::string
        buffer.size(),
        type,
        src,
        comm_,
        req.create());
    comm_mutex_t::instance().unlock();
    if (ret) exit(ret);
    return req;
  }

};

} // namespace

#endif
