#ifndef BLOCK_CONNECT_HPP
#define BLOCK_CONNECT_HPP

#include "config.hpp"
#include "block_ijk.hpp"

#include <map>
#include <vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Block class for computing block connectivity
////////////////////////////////////////////////////////////////////////////////
class mesh_ijk_t {

  std::vector<real_t> vertices_;
  std::vector<block_ijk_t> blocks_;
  
  std::map<std::vector<int_t>, std::vector<int_t>> block_connectivity_;
  
  int_t num_dims_ = -1;

public:

  mesh_ijk_t(int_t num_dims) : num_dims_(num_dims) {}

  mesh_ijk_t(const mesh_ijk_t &) = delete;

  int_t num_roots() const { return blocks_.size(); }
  int_t num_vertices() const { return vertices_.size() / num_dims_; }

  void reserve_blocks(int_t nblock) { blocks_.reserve(nblock); }
  void reserve_vertices(int_t nverts) { vertices_.reserve(num_dims_ * nverts); }
  
  void resize_vertices(int_t nverts) { vertices_.resize(num_dims_ * nverts); }

  void clear_blocks() {
    block_connectivity_.clear();
    blocks_.clear();
  }

  int_t add_vertex(const real_t * vert);
  void set_vertex(int_t id, const real_t * vert);

  const real_t * vertex(int_t id) const
  { return &vertices_[id*num_dims_]; }

  block_ijk_t & add_block(
      const long_t * dims,
      const int_t * verts,
      size_t nverts);

  void compute_neighbors();

  bool partition(int_t nparts, int_t min_size = 1);

  const block_ijk_t & block(int_t b) const { return blocks_[b]; }
  block_ijk_t & block(int_t b) { return blocks_[b]; }

  const real_t * coordinates() const { return &vertices_[0]; }

  template<typename P>
  void traverse(P && p) {
    for (auto & b : blocks_)
      b.descend(std::forward<P>(p));
  }
  
  size_t count_leaves() const;
  
  void get_leaves(std::vector<const block_ijk_t*> & list) const;
  std::vector<const block_ijk_t*> get_leaves() const;

  void number();

  const block_ijk_t * find_largest_block() const;
};

} // namespace

#endif // BLOCK_CONNECT_HPP
