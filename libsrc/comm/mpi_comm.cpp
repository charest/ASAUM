#include "mpi_comm.hpp"
#include "utils/errors.hpp"

#include <iostream>
#include <vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Requests
////////////////////////////////////////////////////////////////////////////////

/// create one
MPI_Request * mpi_request_t::create()
{ 
  requests_.resize(requests_.size() + 1, MPI_REQUEST_NULL);
  return &requests_.back();
}
  
/// transfer one
void mpi_request_t::transfer(mpi_request_t & req)
{
  auto & reqs = req.requests_;
  requests_.insert(requests_.end(), reqs.begin(), reqs.end());
  reqs.clear();
}

/// wait on a request
int mpi_request_t::wait() {
  if (!requests_.empty()) {
    comm_mutex_t::instance().lock();
    auto res = MPI_Waitall(requests_.size(), requests_.data(), MPI_STATUS_IGNORE);
    comm_mutex_t::instance().unlock();
    requests_.clear();
    return res;
  }
  return MPI_SUCCESS;
}

/// check if a request is finished
bool mpi_request_t::finished() {
  if (!requests_.empty()) {
    int flag;
    comm_mutex_t::instance().lock();
    auto res = MPI_Testall(
        requests_.size(),
        requests_.data(),
        &flag,
        MPI_STATUSES_IGNORE);
    if (res) THROW_ERROR("MPI_Testall returned failure");
    return flag;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Function
////////////////////////////////////////////////////////////////////////////////

/// create one
MPI_Op * mpi_function_t::create()
{
  op_ = std::shared_ptr<MPI_Op>(
      new MPI_Op,
      [](MPI_Op * op){ 
        MPI_Op_free(op);
        delete op;
      } );
  return op_.get();
}


////////////////////////////////////////////////////////////////////////////////
// default constructor
////////////////////////////////////////////////////////////////////////////////
mpi_comm_t::mpi_comm_t(int argc, char * argv[], bool is_threaded) :
  comm_(MPI_COMM_WORLD), top_(true) 
{
  comm_mutex_t::instance().lock();

  if (is_threaded) {
    int prov;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov);
    if (prov != MPI_THREAD_MULTIPLE) {
      MPI_Comm_rank(comm_, &rank_);
      if (rank_ == 0) {
        std::cout << "Asking for  MPI_THREAD_MULTIPLE, but it is not available!" << std::endl;
        std::cout << "Recompile MPI with multi-threaded support to enable." << std::endl;
      }
      exit(consts::failure);
    }
  }
  else {
    MPI_Init(&argc, &argv);
  }

  MPI_Comm_rank(comm_, &rank_);
  MPI_Comm_size(comm_, &size_);
  MPI_Comm_group(comm_, &group_);
  
  int len;
  MPI_Get_processor_name(hostname_, &len);
  
  comm_mutex_t::instance().unlock();
}

////////////////////////////////////////////////////////////////////////////////
// contruct subcommmunicator as a duplicate
////////////////////////////////////////////////////////////////////////////////
mpi_comm_t::mpi_comm_t(const mpi_comm_t & comm)
{
  comm_mutex_t::instance().lock();

  MPI_Group parent_group = comm.group_;
  MPI_Comm parent_comm = comm.comm();
    
  if (parent_comm != MPI_COMM_NULL && comm.rank_ >= 0) {

    MPI_Group_excl(
      parent_group,
      0,
      nullptr,
      &group_);
    MPI_Comm_create_group(parent_comm, group_, 0, &comm_);
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
  }
  else {
    rank_ = -1;
    size_ = 0;
  }
  
  int len;
  MPI_Get_processor_name(hostname_, &len);

  comm_mutex_t::instance().unlock();
}


////////////////////////////////////////////////////////////////////////////////
// contruct subcommmunicator from flag
////////////////////////////////////////////////////////////////////////////////
mpi_comm_t::mpi_comm_t(
    const mpi_comm_t & comm,
    bool flag)
{
  
  comm_mutex_t::instance().lock();

  int parent_size = comm.size();
  MPI_Comm parent_comm = comm.comm();
  MPI_Group parent_group = comm.group_;

  std::vector<char> rank_flags(parent_size);
  
  MPI_Allgather(&flag, 1, MPI_CHAR, rank_flags.data(), 1, MPI_CHAR, parent_comm);
  
  int num_ranks = 0;
  for (auto i : rank_flags) if (i) num_ranks++;

  std::vector<int> included_ranks;
  included_ranks.reserve(num_ranks);
  for (int i=0; i<parent_size; ++i)
    if (rank_flags[i])
      included_ranks.emplace_back(i);

  MPI_Group_incl(
    parent_group,
    num_ranks,
    included_ranks.data(),
    &group_);
  
  MPI_Comm_create_group(parent_comm, group_, 0, &comm_);

  if (flag) {
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
  }
  else {
    rank_ = -1;
    size_ = 0;
  }
  
  int len;
  MPI_Get_processor_name(hostname_, &len);

  comm_mutex_t::instance().unlock();
}

////////////////////////////////////////////////////////////////////////////////
// contruct subcommmunicator from flag
////////////////////////////////////////////////////////////////////////////////
mpi_comm_t::mpi_comm_t(
    const mpi_comm_t & comm,
    const std::vector<char> & rank_flags)
{
  comm_mutex_t::instance().lock();

  int parent_size = comm.size();
  MPI_Comm parent_comm = comm.comm();
  MPI_Group parent_group = comm.group_;

  int num_ranks = 0;
  for (auto i : rank_flags) if (i) num_ranks++;

  std::vector<int> included_ranks;
  included_ranks.reserve(num_ranks);
  for (int i=0; i<parent_size; ++i)
    if (rank_flags[i])
      included_ranks.emplace_back(i);

  MPI_Group_incl(
    parent_group,
    num_ranks,
    included_ranks.data(),
    &group_);

  MPI_Comm_create_group(parent_comm, group_, 0, &comm_);
  
  if (rank_flags[comm.rank()]) {
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
  }
  else {
    rank_ = -1;
    size_ = 0;
  }
  
  int len;
  MPI_Get_processor_name(hostname_, &len);

  comm_mutex_t::instance().unlock();
}

////////////////////////////////////////////////////////////////////////////////
// destructor
////////////////////////////////////////////////////////////////////////////////
mpi_comm_t::~mpi_comm_t()
{
  if (comm_ != MPI_COMM_NULL && comm_ != MPI_COMM_WORLD) {
    comm_mutex_t::instance().lock();
    MPI_Comm_free(&comm_);
    comm_mutex_t::instance().unlock();
  }
  if (group_ != MPI_GROUP_NULL){
    comm_mutex_t::instance().lock();
    MPI_Group_free(&group_);
    comm_mutex_t::instance().unlock();
  }
  if (top_) shutdown();
}
  
////////////////////////////////////////////////////////////////////////////////
// shutdown
////////////////////////////////////////////////////////////////////////////////
void mpi_comm_t::shutdown() { 
  mpi_type_registry_t::instance().data_types.clear();
  comm_mutex_t::instance().lock();
  MPI_Finalize();
  comm_mutex_t::instance().unlock();
}
  
void mpi_comm_t::exit(int i) {
  shutdown();
  ::exit(i);
}

////////////////////////////////////////////////////////////////////////////////
// create an op
////////////////////////////////////////////////////////////////////////////////
mpi_function_t mpi_comm_t::create_op(
    mpi_function_t::fun_type fun,
    bool commutes )
{
  mpi_function_t mpi_fun; 
  comm_mutex_t::instance().lock();
  auto ret = MPI_Op_create(fun, commutes, mpi_fun.create());
  comm_mutex_t::instance().unlock();
  if (ret) exit(ret);
  return mpi_fun;
}

////////////////////////////////////////////////////////////////////////////////
// create a barrier
////////////////////////////////////////////////////////////////////////////////
void mpi_comm_t::barrier() 
{
  comm_mutex_t::instance().lock();
  MPI_Barrier(comm_);
  comm_mutex_t::instance().unlock();
}

////////////////////////////////////////////////////////////////////////////////
// probe for a message
////////////////////////////////////////////////////////////////////////////////
bool mpi_comm_t::probe(int source, int tag)
{
  MPI_Status status;
  int flag = 0;
  comm_mutex_t::instance().lock();
  auto ret = MPI_Iprobe(source, tag, comm_, &flag, &status);
  comm_mutex_t::instance().unlock();
  if (ret) exit(ret);
  return flag;
}

bool mpi_comm_t::probe()
{
  MPI_Status status;
  int flag = 0;
  comm_mutex_t::instance().lock();
  auto ret = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &flag, &status);
  comm_mutex_t::instance().unlock();
  if (ret) exit(ret);
  return flag;
}

bool mpi_comm_t::probe_any(int & rank, int & tag)
{
  MPI_Status status;
  int flag = 0;
  comm_mutex_t::instance().lock();
  auto ret = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &flag, &status);
  comm_mutex_t::instance().unlock();
  if (ret) exit(ret);
  rank = status.MPI_SOURCE;
  tag = status.MPI_TAG;
  return flag;
}

bool mpi_comm_t::probe_tag(int & rank, int tag)
{
  MPI_Status status;
  int flag = 0;
  comm_mutex_t::instance().lock();
  auto ret = MPI_Iprobe(MPI_ANY_SOURCE, tag, comm_, &flag, &status);
  comm_mutex_t::instance().unlock();
  if (ret) exit(ret);
  rank = status.MPI_SOURCE;
  return flag;
}


} // namespace
