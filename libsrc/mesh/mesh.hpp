#ifndef MESH_HPP
#define MESH_HPP

#include "config.hpp"
#include "mesh_block.hpp"
#include "mesh_node.hpp"
#include "mesh_boundary.hpp"

#include "comm/comm_map.hpp"
#include "comm/mpi_comm.hpp"

#include <vector>

namespace prl {


struct amr_flags_t;
struct comm_map_block_t;
struct face_id_t;
class lua_result_t;
class mesh_struct_block_t;
class mesh_unstruct_block_t;
class mesh_struct_node_t;
class mesh_unstruct_node_t;

/// The different distribution types
enum class distribution_t {
  sequential,
  cyclic,
  hostname
};

////////////////////////////////////////////////////////////////////////////////
/// Base mesh class
////////////////////////////////////////////////////////////////////////////////
class mesh_t {

protected:
  
  int_t num_dims_ = 0;

  int_t max_refinement_level_ = 5;
  
  std::vector<int_t> partitioning_;

  std::vector<int_t> block_offsets_;
  
  mpi_comm_t & comm_;
  
  std::vector<std::unique_ptr<mesh_node_t>> roots_;
  std::vector<mesh_block_t*> blocks_;

  std::vector<std::unique_ptr<mesh_boundary_t>> boundaries_;
  std::map<std::string, int_t> boundary_map_;
  
  std::unique_ptr<comm_map_t> comm_maps_; // Must come before any objects that
                                          // use this map with a comm queue.
                                          // Since the comm queue cannot destroy
                                          // itself without the map.  An
                                          // alternative would be to have a
                                          // user defined destructor.

public:

  mesh_t(
      lua_result_t input,
      mpi_comm_t & comm,
      bool is_structured,
      const std::vector<int_t> & parts);

  virtual bool is_structured() const = 0;
  
  int_t num_dims() const { return num_dims_; }
  
  int_t num_blocks() const { return blocks_.size(); }

  int_t tot_blocks() const
  { return block_offsets_.empty() ? 0 : block_offsets_.back(); }

  bool owns_block(int_t bid) const;
  int_t block_owner(int_t bid) const;

  const auto & block_offsets() const { return block_offsets_; }

  const auto & boundary_map() const { return boundary_map_; }

  int_t block_global2local(int_t bid)
  { return bid - block_offsets_[comm_.rank()]; }

  std::vector<int_t> block_counts() const;
  std::vector<int_t> block_ids() const;
  
  mesh_block_t * block(int_t i)
  {
    if (i<0 || i>=num_blocks()) {
      THROW_ERROR(
          "Asking for block " << i << " on processor "
          << comm_.rank() << " but only " << num_blocks() << " exist.");
    }
    return blocks_[i];
  }
  
  const mesh_block_t * block(int_t i) const
  {
    if (i<0 || i>=num_blocks()) {
      THROW_ERROR(
          "Asking for block " << i << " on processor "
          << comm_.rank() << " but only " << num_blocks() << " exist.");
    }
    return blocks_[i];
  }
  
  mesh_struct_block_t * struct_block(int_t i);
  mesh_unstruct_block_t * unstruct_block(int_t i);
  const mesh_struct_block_t * struct_block(int_t i) const;
  const mesh_unstruct_block_t * unstruct_block(int_t i) const;

  int_t num_cells() const {
    int_t num = 0;
    for (const auto & b : blocks_) num += b->num_owned_cells();
    return num;
  }

  auto comm_map() { return comm_maps_.get(); }
  
  void set_max_refinement_level(int_t i) { max_refinement_level_ = i; }

  mesh_node_t * root(int_t i) { return roots_[i].get(); }
  void reserve_roots(size_t n)
  { roots_.reserve(n); }

  mesh_struct_node_t * struct_root(int_t i);
  mesh_unstruct_node_t * unstruct_root(int_t i);
  const mesh_struct_node_t * struct_root(int_t i) const;
  const mesh_unstruct_node_t * unstruct_root(int_t i) const;

  void add_root(std::unique_ptr<mesh_node_t> node)
  { roots_.emplace_back( std::move(node) ); }

  void build_geometry(bool faces, bool cells);
  void exchange_geometry();
  
  void build_halo(int_t num_ghost=1, bool with_corners=false);
  virtual void _build_halo(
      int_t num_ghost,
      bool with_corners,
      std::vector<comm_map_block_t> & comm_maps) = 0;

  void build_boundaries();

  void prune_connectivity()
  { for (const auto & b : blocks_) b->prune_connectivity(); }

  virtual void number() = 0;
  
#ifdef HAVE_EXODUS
  virtual void output(exo_writer_t & exo, int_t blk, const std::string & lbl)
  { blocks_[blk]->output(exo, lbl); }
#endif

  /// Find which cell a position is located in
  bool find_cell(const real_t * x, int_t & block_id, int_t & cell_id);
  
  bool find_cell(long_t gid, int_t & block_id, int_t & cell_id);

  void ghost_exchange();

  void number_disjoint(int_t dim);
  
  size_t count_leaves() const;
  void number_leaves();
  
  void get_leaves(std::vector<mesh_node_t*> & list) const;
  std::vector<mesh_node_t*> get_leaves() const;
  
  void get_leaf_blocks(std::vector<mesh_block_t*> & list) const;
  std::vector<mesh_block_t*> get_leaf_blocks() const;
  
  void fixup(
      std::vector<mesh_node_t*> & blocks,
      std::vector<amr_flags_t> & flags);
  
  bool check(
      std::vector<mesh_node_t*> & blocks,
      std::vector<amr_flags_t> & flags,
      bool first_pass);

  bool inner_check(
      std::vector<mesh_node_t*> & blocks,
      std::vector<amr_flags_t> & flags,
      bool first_pass);
  
  bool final_check(
      std::vector<mesh_node_t*> & blocks,
      std::vector<amr_flags_t> & flags);

  void change_connect(
      std::vector<mesh_node_t*> & blocks,
      const std::vector<amr_flags_t> & flags);

  bool refine_tree(
      std::vector<mesh_node_t*> & blocks,
      const std::vector<amr_flags_t> & flags,
      amr_callback_t callback);
  
  bool coarsen_tree(
      std::vector<mesh_node_t*> & blocks,
      const std::vector<amr_flags_t> & flags,
      amr_callback_t callback);

  void post_amr();

  bool amr(
      const std::vector<amr_flags_t> & flags,
      amr_callback_t rcb = {},
      amr_callback_t ccb = {});

  void refine(const bool * flags);
  void refine(int_t * levels);

  void initial_refinement(const std::vector<int_t> & levels);

  void find_neighbors();

  virtual ~mesh_t() {}
};

////////////////////////////////////////////////////////////////////////////////
/// utility function
////////////////////////////////////////////////////////////////////////////////
std::unique_ptr<mesh_t> make_mesh(
    lua_result_t inputs,
    mpi_comm_t & mycomm,
    const std::vector<int_t> & parts);


} // namespace

#endif
