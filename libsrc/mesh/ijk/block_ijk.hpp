#ifndef BLOCK_HPP
#define BLOCK_HPP

#include "config.hpp"
#include "block_ijk_neighbor.hpp"
#include "block_ijk_vertex.hpp"
#include "range.hpp"
#include "sector.hpp"
#include "math/subdivide.hpp"
#include "utils/result.hpp"

#include <functional>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>

namespace prl {

class mesh_ijk_t;
class transform_t;

////////////////////////////////////////////////////////////////////////////////
/// Block class for computing block connectivity
////////////////////////////////////////////////////////////////////////////////
class block_ijk_t {

public:

  struct compare_t {
    bool operator()(const std::vector<int_t> & a, const std::vector<int_t> & b) const
    {
      if (a.size() == b.size())
        return a < b;
      else
        return a.size() > b.size();
    }
  };
  
  using child_list_t = std::vector<std::unique_ptr<block_ijk_t>>;

  using int3_t = std::array<int_t, 3>;
  using long3_t = std::array<long_t,3>;
  
private:

  static sector_map_t<int3_t> sector_ctm_[3];

  int_t id_ = -1;
  int_t label_ = -1;
  int_t split_ = -1;
  int3_t level_ = {0, 0, 0};
  int_t num_dims_ = 0;

  range_ijk_t range_;
  long3_t dims_ = {1, 1, 1};

  block_ijk_t * parent_ = nullptr;
  child_list_t children_;

  std::vector<int_t> vertex_ids_;
  
  sector_map_t<std::vector<int_t>> sector2verts_;
  std::map<std::vector<int_t>, sector_t, compare_t> verts2sector_;
  
  sector_map_t<block_ijk_vertex_t> vertices_;
  std::map<int_t, sector_t> vertex_map_;

  sector_map_t<std::vector<block_ijk_neighbor_t>> neighbors_;
  std::vector<block_ijk_shared_t> shared_;
  
  std::function<void(const real_t *, const long_t*, const long_t *, real_t*)> spacing_;

  mesh_ijk_t & mesh_;

public:
  
  block_ijk_t(
      int_t id,
      const long_t * end, 
      const int_t * verts,
      size_t ndims,
      size_t nverts,
      mesh_ijk_t & mesh,
      const long_t * begin = nullptr);
  
  block_ijk_t(
      block_ijk_t * parent,
      const range_ijk_t & range, 
      const int_t * verts,
      size_t ndims,
      size_t nverts,
      mesh_ijk_t & mesh);

  bool operator==(const block_ijk_t & other) const
  { return this == &other; }

  size_t size() const
  { return std::accumulate(&dims_[0], &dims_[num_dims_], 1, std::multiplies<size_t>{}); }

  int_t split_dim() const { return split_; }
  int_t id() const { return id_; }
  int_t label() const { return label_; }
  
  const auto & levels() const { return level_; }
  int_t level(int_t dim) const { return level_[dim]; }

  bool is_active() const
  { return id_>-1; }

  bool is_root() const
  { return !parent_; }

  const auto begin() const { return &range_.begin[0]; }
  const auto begin(int_t i) const { return range_.begin[i]; }

  const auto end() const { return &range_.end[0]; }
  const auto end(int_t i) const { return range_.end[i]; }

  int_t num_vertices() const { return vertex_ids_.size(); }

  const block_ijk_t * parent() const { return parent_; }
  int_t parent_id() const { return parent_ ? parent_->id() : id_; }

  const block_ijk_t * root() const { return parent_ ? parent_->root() : this; }
  
  const long_t * dims() const { return &dims_[0]; }
  long_t dim(int_t i) const { return dims_[i]; }
  
  void connect(int_t i, int_t j, int_t k, std::vector<int_t> conn);
  
  void add_vertex(sector_t sec, block_ijk_vertex_t v);

  void orient(
      const std::vector<int_t> & sorted_vs,
      block_ijk_t & donor);

  int_t num_dims() const { return num_dims_; }
  
  int_t num_cells() const 
  { return std::accumulate(dims_.begin(), dims_.end(), 1, std::multiplies<int_t>()); }

  int_t num_neighbors() const 
  { 
    int_t n=0;
    for (const auto & neighs : neighbors_)
      n += neighs.second.size();
    return n;
  }

  void make_relative(range_ijk_t & rng) const;

  template<typename T>
  void make_relative(T * pos) const
  {
    for (int_t dim=0; dim<num_dims_; ++dim) {
      pos[dim] -= range_.begin[dim];
    }
  }
  
  void add_internal_neighbor(
      block_ijk_t & neigh,
      const sector_t & my_sector,
      range_ijk_t my_range);
  
  void add_external_neighbor(
      block_ijk_t & neigh,
      sector_t my_sector,
      range_ijk_t my_range,
      transform_t me2donor_cells,
      transform_t me2donor_verts);

  template<typename...ARGS>
  void add_neighbor(const sector_t & sector, ARGS &&... args)
  { neighbors_[sector].emplace_back(std::forward<ARGS>(args)...); }
  
  template<typename...ARGS>
  void add_shared(ARGS &&... args)
  { shared_.emplace_back(std::forward<ARGS>(args)...); }


  bool has_neighbor_at(int_t i, int_t j, int_t k) const
  { return neighbors_.count({i, j, k}); }

  bool has_neighbor_at(sector_t sec) const
  { return neighbors_.count(sec); }
  
  bool has_neighbor_at(bool is_hi, int_t dim) const
  { return neighbors_.count( {is_hi, dim} ); }
  
  bool has_neighbor_at(sector_t sec, const block_ijk_t & n) const;
  
  const auto & neighbors_at(int_t i, int_t j, int_t k) const
  { return neighbors_.at({i, j, k}); }

  const auto & neighbors_at(sector_t sec) const
  { return neighbors_.at(sec); }

  const auto & neighbors_at(bool is_hi, int_t dim) const
  { return neighbors_.at( {is_hi, dim} ); }

  const auto & neighbors() const { return neighbors_; }
  const auto & shared() const { return shared_; }

  const std::vector<int_t> & vertices() const { return vertex_ids_; }

  std::vector<real_t> coords() const;
  void coords(real_t * xs) const;
  
  void set_linear_spacing();

  result_t<sector_t> sector(std::vector<int_t> vs) const;

  void spacing(const real_t * block_coords, const long_t* is, const long_t* dims, real_t* xs) const
  { spacing_(block_coords, is, dims, xs); }

  const auto & verts2sector() const { return verts2sector_; }

  bool split_block(size_t min);

  void split_block(const range_ijk_t * ranges, int_t nranges);

  const range_ijk_t & range() const { return range_; }
  range_ijk_t range(const sector_t & sector) const;
  range_ijk_t range_wrt_root(const sector_t & sector) const;

  void number(int_t & counter);
  
  void get_leaves(std::vector<const block_ijk_t*> & list) const;
  std::vector<const block_ijk_t*> get_leaves() const;
  
  result_t<child_list_t::const_iterator> find_child(const block_ijk_t * c) const;

  const block_ijk_t * left(const block_ijk_t * c) const;
  
  const block_ijk_t * right(const block_ijk_t * c) const;
 
  template<typename P>
  void descend(P && p) const
  {
    if (children_.empty()) {
      std::forward<P>(p)(*this);
    }
    else {
      for (auto & c : children_)
        c->descend(std::forward<P>(p));
    }
  }

  size_t count_leaves() const;
  
  const block_ijk_t * find_largest_block() const;
  
  friend std::ostream &operator<<( std::ostream &out, const block_ijk_t & n ) {
    auto nd = n.num_dims();
    out << " id=" << n.id() << ", level=(";
    for (int_t d=0; d<nd-1; ++d) out <<  n.level(d) << ", ";
    if (nd>0) out <<  n.level(nd-1);
    out << "), range=[(";
    for (int_t d=0; d<nd-1; ++d) out << n.begin(d) << ", ";
    if (nd>0) out << n.begin(nd-1) << "), (";
    for (int_t d=0; d<nd-1; ++d) out << n.end(d) << ", ";
    if (nd>0) out << n.end(nd-1);
    out << ")]";
    return out;
  }
  
private:
  
  void _find_sectors(
      const block_ijk_t * root,
      const sector_t & sector,
      std::vector<const block_ijk_t *> & list) const;

  void _find_neighbors(
      const block_ijk_t * root,
      const sector_t & sector,
      std::vector<const block_ijk_t *> & list,
      std::array<bool,3> & mask) const;
  
public:
  
  void find_sectors(
      const sector_t & sector,
      std::vector<const block_ijk_t *> & list) const;
  std::vector<const block_ijk_t *> find_sectors(const sector_t & sector) const;
  
  void find_neighbors(
      const sector_t & sector,
      std::vector<const block_ijk_t *> & list) const;
  std::vector<const block_ijk_t *> find_neighbors(const sector_t & sector) const;
  
  void find_neighbors_descend(
      const range_ijk_t & range,
      const sector_t & sector,
      std::vector<const block_ijk_t *> & list) const;
  
  std::vector<const block_ijk_t *> find_neighbors_descend(
      const range_ijk_t & range,
      const sector_t & sector) const;

  void find_neighbors_descend(
      const sector_t & sector,
      std::vector<const block_ijk_t *> & list) const;
  std::vector<const block_ijk_t *> find_neighbors_descend(const sector_t & sector) const;
  
};

} // namespace

#endif
