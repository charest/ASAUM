#ifndef MESH_STRUCT_NODE_HPP
#define MESH_STRUCT_NODE_HPP

#include "mesh_node.hpp"
#include "ijk/ctm.hpp"
#include "ijk/sector.hpp"
#include "ijk/transform.hpp"
#include "utils/result.hpp"

namespace prl {

class mesh_struct_block_t;
class mesh_struct_node_t;

////////////////////////////////////////////////////////////////////////////////
/// Neighbor info class 
////////////////////////////////////////////////////////////////////////////////
class struct_neighbor_t {
  const mesh_node_t * neigh_;
  ctm_t me2donor_cells_, me2donor_verts_;
  range_ijk_t range_;

public:

  struct_neighbor_t(
      const mesh_node_t * neigh,
      const ctm_t & trans_cells,
      const ctm_t & trans_verts,
      const range_ijk_t & range) :
    neigh_(neigh),
    me2donor_cells_(trans_cells),
    me2donor_verts_(trans_verts),
    range_(range)
  {}

  const mesh_node_t * node() const { return neigh_; }

  int_t id() const { return neigh_->id(); }

  bool is_used() const { return neigh_->is_used(); }
  
  const auto & range() const { return range_; }
  
  transform_t tm_me2donor_cells() const { return me2donor_cells_; }
  transform_t tm_me2donor_verts() const { return me2donor_verts_; }
  
  const ctm_t & ctm_me2donor_cells() const { return me2donor_cells_; }
  const ctm_t & ctm_me2donor_verts() const { return me2donor_verts_; }

  
  void me2donor_cells(int_t * is) const
  { me2donor_cells_.transform(is); }
  
  void me2donor_verts(int_t * is) const
  { me2donor_verts_.transform(is); }
  
  bool vert_is_inside(const int_t * pos) const
  { return range_.is_inside(pos, neigh_->num_dims(), true); }

};

////////////////////////////////////////////////////////////////////////////////
/// Shared info class 
////////////////////////////////////////////////////////////////////////////////
class struct_shared_t {
  const mesh_node_t * neigh_;
  sector_t sector_;
  range_ijk_t range_;

public:

  struct_shared_t(
      const mesh_node_t * neigh,
      const sector_t & sector,
      const range_ijk_t & my_range) : 
    neigh_(neigh), sector_(sector), range_(my_range)
  {}

  const mesh_node_t * node() const { return neigh_; }

  const auto & sector() const { return sector_; }
  const auto & range() const { return range_; }
  
  int_t id() const { return neigh_->id(); }
  bool is_used() const { return neigh_->is_used(); }
};


////////////////////////////////////////////////////////////////////////////////
/// The main  messh node 
////////////////////////////////////////////////////////////////////////////////
class mesh_struct_node_t : public mesh_node_t {
public:

  using neighbor_pair_t = 
    std::pair<const mesh_struct_node_t *, const struct_neighbor_t *>;

private:
  
  int_t split_ = -1;
  
  std::array<long_t,3> dims_ = {1, 1, 1};
  range_ijk_t range_;
  
  sector_map_t<std::vector<struct_neighbor_t>> root_neighbors_;

  sector_map_t<std::vector<struct_neighbor_t>> neighbors_;
  std::vector<struct_shared_t> shared_;

 
public:

  mesh_struct_node_t(const long_t * dims, int_t ndims, int_t id);
  
  mesh_struct_node_t(mesh_struct_node_t * parent, int_t split, int_t sector, long_t midp);

  const long_t * dims() const { return &dims_[0]; }
  long_t dim(int_t i) const { return dims_[i]; }
  
  const auto begin() const { return &range_.begin[0]; }
  const auto begin(int_t i) const { return range_.begin[i]; }

  const auto end() const { return &range_.end[0]; }
  const auto end(int_t i) const { return range_.end[i]; }

  int_t split_dim() const { return split_; }

  mesh_struct_block_t * struct_block();

  mesh_struct_node_t * sibling() const;
  mesh_struct_node_t * left(const mesh_struct_node_t * c) const;
  mesh_struct_node_t * right(const mesh_struct_node_t * c) const;

  template<typename...Args>
  void add_neighbor(const sector_t & sector, Args&&...args)
  { neighbors_[sector].emplace_back(std::forward<Args>(args)...); }
  
  template<typename...Args>
  void add_root_neighbor(const sector_t & sector, Args&&...args)
  { root_neighbors_[sector].emplace_back(std::forward<Args>(args)...); }

  void add_internal_neighbor(
      mesh_struct_node_t * neigh,
      const sector_t & sector);

  void add_external_neighbor(
    mesh_struct_node_t * neigh,
    sector_t my_sector,
    transform_t me2neigh_cells,
    transform_t me2neigh_verts);
  
  template<typename...Args>
  void add_shared(Args&&...args)
  { shared_.emplace_back(std::forward<Args>(args)...); }
  
  bool has_neighbor_at(sector_t sec) const
  { return neighbors_.count(sec); }
  
  const auto & neighbors_at(sector_t sec) const
  { return neighbors_.at(sec); }

  const auto & neighbors() const { return neighbors_; }
  const auto & shared() const { return shared_; }
  
  const range_ijk_t & range() const { return range_; }
  range_ijk_t range(const sector_t & sector) const;
  range_ijk_t range_wrt_root(const sector_t & sector) const;

  void make_relative(long_t * pos) const;
  void make_relative(range_ijk_t & rng) const;
  
  bool first_sibling_check(amr_flags_t *) override;
  
  bool ratio_check(amr_flags_t *, bool) override;

  //bool final_check_1(amr_flags_t *, const amr_flags_t *) override;

  //bool final_check_2(amr_flags_t *, const amr_flags_t *) override;
 
  //void change_connect(const amr_flags_t *) override;
  
  void insert_split(int_t dim);
  virtual void merge();
  
  bool refine_tree(const amr_flags_t & flags) override;
  bool coarsen_tree(const amr_flags_t & flags) override;
  
  void refine_block(const amr_flags_t & flags) override;
  void coarsen_block() override;

  void clear_neighbors() override;
  void find_neighbors() override;
  
  void find_neighbors(
      const sector_t & sector,
      std::vector<neighbor_pair_t> & list) const;

  void find_neighbors(
      const mesh_struct_node_t * start,
      const sector_t & sector,
      std::array<bool,3> & mask,
      std::vector<neighbor_pair_t> & list) const;
  
  void find_neighbors_descend(
      const int_t * level,
      const range_ijk_t & range,
      const sector_t & sector,
      std::vector<neighbor_pair_t> & list) const;
  
  void find_siblings(
      int_t dim,
      std::vector<const mesh_struct_node_t*> & list) const;
  void find_siblings_descend(
      int_t dim,
      bool go_right,
      std::vector<const mesh_struct_node_t*> & list) const;

};


} // namespace

#endif
