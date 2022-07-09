#include "config.hpp"
#include "cast.hpp"
#include "trace.hpp"

#include "comm/mpi_comm.hpp"

#include <cassert>
#include <fstream>
#include <numeric>

namespace prl {

/// This is the first instance used, to start the top clock
auto & obj = tracer_t::instance();
  
////////////////////////////////////////////////////////////////////////////////
/// truncate 
////////////////////////////////////////////////////////////////////////////////
std::string gathered_trace_node_t::pretty_name(size_t len) const
{
  if (name_.size() > len) {
    auto name = name_;
    name.resize(len);
    return name;
  }
  else {
    return name_;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// De-serialize trace data
////////////////////////////////////////////////////////////////////////////////
void gathered_trace_node_t::unpack(
    const byte_t * & buffer,
    int_t rank)
{
  size_t name_size;
  uncast( buffer, 1, &name_size );

  name_.clear();
  name_.resize(name_size);
  uncast( buffer, name_size, name_.data() ); 
  
  uncast( buffer, 1, &elapsed_[rank] );
  uncast( buffer, 1, &counter_[rank] );

  size_t num_child;
  uncast( buffer, 1, &num_child );

  for (size_t c=0; c<num_child; ++c) {
    
    auto child_buffer = buffer;
    uncast( child_buffer, 1, &name_size );

    std::string child_name;
    child_name.resize(name_size);
    uncast( child_buffer, name_size, child_name.data() );

    auto it = children_.find(child_name);
    if (it == children_.end()) {
      auto child = std::make_unique<gathered_trace_node_t>();
      child->unpack(buffer, rank);
      children_.emplace(child->name_, std::move(child));
    }
    else {
      auto & child = it->second;
      child->unpack(buffer, rank);
    }

  } // child
}

////////////////////////////////////////////////////////////////////////////////
/// Get total elapsed
////////////////////////////////////////////////////////////////////////////////
double gathered_trace_node_t::elapsed() const {
  double elapse = 0;
  for (const auto & e : elapsed_)
    elapse += e.second;
  return elapse;
}

////////////////////////////////////////////////////////////////////////////////
/// Get the totall counts
////////////////////////////////////////////////////////////////////////////////
size_t gathered_trace_node_t::counts() const {
  size_t sum = 0;
  for (const auto & e : counter_)
    sum += e.second;
  return sum;
}
  
////////////////////////////////////////////////////////////////////////////////
/// Dump info for this node to screen
////////////////////////////////////////////////////////////////////////////////
void gathered_trace_node_t::print(size_t level, double tot) const {
  auto space = level*2;
  auto elaps = elapsed();
  auto cnts = counts();
  auto len = 36-space;
  auto ss = std::cout.precision();
  std::cout << std::string(level*2, ' ') << std::setw(len) << std::left << pretty_name(len);
  std::cout.setf( std::ios::scientific );
  std::cout << std::setw(12) << std::right << std::setprecision(4) << elaps;
  std::cout.unsetf( std::ios::scientific );
  std::cout.setf( std::ios::fixed );
  std::cout << std::setw(10) << std::right << std::setprecision(2) << elaps / tot * 100.;
  std::cout.unsetf( std::ios::fixed );
  std::cout << std::setw(12) << std::right << cnts << std::endl;
  std::cout.precision(ss);
  for (const auto & child : children_)
    child.second->print(level+1, tot);
}

////////////////////////////////////////////////////////////////////////////////
/// Dump info for this node to file
////////////////////////////////////////////////////////////////////////////////
void gathered_trace_node_t::dump(
    std::ostream & out,
    size_t level,
    double tot,
    size_t nranks) const
{
  out << std::string(level*2, ' ') << name_ << ", ";

  auto ss = out.precision();

  out.setf( std::ios::scientific );

  double tot_elapsed = 0, max_elapsed = 0, min_elapsed = consts::real_max;

  for (size_t r=0; r<nranks; ++r) {
    double elaps = elapsed_.count(r) ? elapsed_.at(r) : 0;
    out << std::setprecision(4) << elaps << ", ";
    tot_elapsed += elaps;
    max_elapsed = std::max( max_elapsed, elaps );
    min_elapsed = std::min( min_elapsed, elaps );
  }
    
  out << std::setprecision(4) << tot_elapsed << ", ";
  out << std::setprecision(4) << tot_elapsed/nranks << ", ";
  out << std::setprecision(4) << max_elapsed << ", ";
  out << std::setprecision(4) << min_elapsed << ", ";
    
  out.unsetf( std::ios::scientific );

  out.setf( std::ios::fixed );
  out << std::setprecision(2) << tot_elapsed / tot * 100.;
  out.unsetf( std::ios::fixed );
    
  for (size_t r=0; r<nranks; ++r) {
    size_t cnt = counter_.count(r) ? counter_.at(r) : 0;
    out << std::setprecision(4) << cnt << ", ";
  }

  out << std::endl;

  out.precision(ss);

  for (const auto & child : children_)
    child.second->dump(out, level+1, tot, nranks);
}

////////////////////////////////////////////////////////////////////////////////
/// truncate 
////////////////////////////////////////////////////////////////////////////////
std::string trace_node_t::pretty_name(size_t len) const
{
  if (name_.size() > len) {
    auto name = name_;
    name.resize(len);
    return name;
  }
  else {
    return name_;
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Dump info for this node to screen
////////////////////////////////////////////////////////////////////////////////
void trace_node_t::print(size_t level, double tot) const {
  auto space = level*2;
  auto len = 36-space;
  auto elaps = elapsed();
  auto cnts = counter_;
  auto ss = std::cout.precision();
  std::cout << std::string(level*2, ' ') << std::setw(len) << std::left << pretty_name(len);
  std::cout.setf( std::ios::scientific );
  std::cout << std::setw(12) << std::right << std::setprecision(4) << elaps;
  std::cout.unsetf( std::ios::scientific );
  std::cout.setf( std::ios::fixed );
  std::cout << std::setw(10) << std::right << std::setprecision(2) << elaps / tot * 100.;
  std::cout.unsetf( std::ios::fixed );
  std::cout << std::setw(12) << std::right << cnts;
  std::cout << std::endl;
  std::cout.precision(ss);
  for (const auto & child : children_)
    child.second->print(level+1, tot);
}

////////////////////////////////////////////////////////////////////////////////
/// Count the size of the tree
////////////////////////////////////////////////////////////////////////////////
void trace_node_t::count(size_t & sz) const {
  sz += sizeof(size_t); // string len
  sz += sizeof(char) * name_.size(); // name
  sz += sizeof(double); // elapsed
  sz += sizeof(size_t); // num children
  sz += sizeof(size_t); // counter
  for (const auto & child : children_)
    child.second->count(sz);
}

////////////////////////////////////////////////////////////////////////////////
/// Pack the tree
////////////////////////////////////////////////////////////////////////////////
void trace_node_t::pack(std::vector<byte_t> & buffer) const
{
  size_t name_size = name_.size();
  cast_insert( &name_size, 1, buffer );
  
  cast_insert(name_.data(), name_.size(), buffer);
  
  auto elaps = elapsed();
  cast_insert(&elaps, 1, buffer);
  cast_insert(&counter_, 1, buffer);
  
  size_t num_child = children_.size();
  cast_insert(&num_child, 1, buffer);

  for (const auto & child : children_)
    child.second->pack(buffer);
}
  
  
////////////////////////////////////////////////////////////////////////////////
/// Dump a programs trace
////////////////////////////////////////////////////////////////////////////////
void tracer_t::dump(mpi_comm_t & comm, bool print, const std::string & filename)
{

  //--- stop the clocks
  if (current_->name() == "main") {
    current_->stop();
  }
  else {
    current_->check("__top__");
  }

  root_.stop();

  //--- pack buffers to send
  auto is_root = comm.is_root();
  auto comm_size = comm.size();
 
  size_t tree_size = 0;
  root_.count(tree_size);

  std::vector<byte_t> sendbuf;
  sendbuf.reserve(tree_size);
  root_.pack(sendbuf);

  if (sendbuf.size() != tree_size)
    THROW_ERROR("Buffer size does not match counts!");

  //--- send buffers to root
  if (is_root) {
    
    std::vector<int> recvcounts(comm_size);
    comm.gather(static_cast<int>(tree_size), recvcounts, 0);

    std::vector<int> recvdispls(comm_size+1);
    recvdispls[0] = 0;
    std::partial_sum(recvcounts.begin(), recvcounts.end(), &recvdispls[1]);

    std::vector<byte_t> recvbuf(recvdispls[comm_size]);
    comm.gatherv(sendbuf, recvbuf, recvcounts, recvdispls, 0);

    gathered_trace_node_t gathered_trace;

    for (int_t r=0; r<comm_size; ++r) {
      auto start = recvdispls[r];
      const auto * buf_start = &recvbuf[start];
      gathered_trace.unpack( buf_start, r );
      assert( buf_start == &recvbuf[recvdispls[r+1]] && "problem unpacking trace" );
    }
  
    // dump trace to screen
    if (print) {
      std::cout << std::string(70, '=') << std::endl;
      std::cout << "BEGIN TRACE" << std::endl;
      std::cout << std::string(36, ' ');
      std::cout << std::setw(12) << std::right << "secs";
      std::cout << std::setw(10) << std::right << "% tot";
      std::cout << std::setw(12) << std::right << "count" << std::endl;
      
      size_t level = 0;
      auto total = gathered_trace.elapsed();
      gathered_trace.print(level, total);

      std::cout << std::endl;
      std::cout << "Note: Times and counts are summed across all ranks." << std::endl;
      std::cout << std::endl;
      std::cout << "END TRACE" << std::endl;
      std::cout << std::string(70, '=') << std::endl;
    } // print

    // dump trace to file
    if (!filename.empty()) {
      std::ofstream file(filename);

      if (file.good()) {
        file << "path, ";
        for (int_t r=0; r<comm_size; ++r)
          file << "time rank " << r << " (secs), ";
        file << " total time (secs), average per rank (secs), ";
        file << " max per rank (secs), min per rank (secs), % total";
        for (int_t r=0; r<comm_size; ++r)
          file << "count rank " << r << ", ";
        file << std::endl;
        size_t level = 0;
        auto total = gathered_trace.elapsed();
        gathered_trace.dump(file, level, total, comm_size);
      }
      else {
        std::cout << "Couldn't create trace file '" << filename << "'." << std::endl;
        comm.exit(consts::failure);
      }

    } // file

  }
  else {
    comm.gather(static_cast<int>(tree_size), 0);
    comm.gatherv(sendbuf, 0);
  }

}

} // namespacce
