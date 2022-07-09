#ifndef STRUCTURED_BLOCK_HPP
#define STRUCTURED_BLOCK_HPP

#include "mesh_block.hpp"
#include "ijk/sector.hpp"
#include "face_data.hpp"
#include "logical/logical_iterators.hpp"
#include "logical/logical_flags.hpp"

#include <memory>
#include <unordered_map>

namespace prl {

class amr_flags_t;
class block_ijk_t;
class mesh_struct_boundary_t;
class mesh_struct_node_t;
class struct_neighbor_t;
class struct_shared_t;
struct entity_info_t;
struct range_ijk_t;
  

////////////////////////////////////////////////////////////////////////////////
/// Structured mesh block
////////////////////////////////////////////////////////////////////////////////
class mesh_struct_block_t : public mesh_block_t {
public:

  int_t region_ = -1;
  std::vector<long_t> dims_;
  logical_block_t iterator_;
  logical_flags_t boundary_flags_;

  struct label_pair_t { 
    bool is_hi=false;
    int_t dim=0; 
    label_pair_t(bool hi, int_t d) : is_hi(hi), dim(d) {}
  };
  
  std::vector<int_t> boundary_side2id_[3][2];
  std::map<int_t, label_pair_t> boundary_id2side_;
  
  std::vector<intrablock_face_t> intrablock_faces_;
  std::vector<int_t> intrablock_cell_to_cell_;
  crs_t<int_t, int_t> intrablock_cell_to_intrablock_faces_;

  std::map<int_t, std::map<std::vector<int_t>, int_t> > sorted_verts2entities_;

  std::map< int_t, std::map<int_t, std::vector<int_t>> > cells2ghost_;
  std::map< int_t, std::map<int_t, int_t> > entities2ghost_;
  std::map< int_t, std::map<int_t, int_t> > ghost2entities_;
  std::map< std::pair<int_t, int_t>, crs_t<int_t, int_t> > ghost_connectivity_;
  
  std::map<int_t, std::vector<long_t>> ghost_local2global_;
  std::map<int_t, std::map<long_t, int_t>> ghost_global2local_;

  mesh_struct_block_t(
      int_t region,
      const long_t * dims,
      int_t ndims);
  mesh_struct_block_t(const block_ijk_t & block);
  
  mesh_struct_node_t * my_node();
  const mesh_struct_node_t * my_node() const;

  int_t count_vertices() const;
  void number_vertices(
      const std::vector<int_t> & block_offsets,
      long_t vert_start,
      std::vector<entity_info_t> & shared_info,
      std::vector<entity_info_t> & ghost_info);
  
  void extract_boundary(int_t boundary_id, std::vector<boundary_face_t> & faces) override;
  
  bool install_boundary(
      const block_ijk_t & block_ijk,
      const mesh_struct_boundary_t & bnd);
  bool install_boundaries(
      const block_ijk_t & block_ijk,
      const std::vector<std::unique_ptr<mesh_boundary_t>> & boundaries);
  
  void add_boundary_mapping(bool is_hi, int_t dir, int_t bid);
  void add_boundary_mappings(bool is_hi, int_t dir, const std::vector<int_t> labels);
  bool find_boundary_from_id(int_t bid, bool & is_hi, int_t & dir) const;

  const std::vector<int_t> & boundary_ids(bool is_hi, int_t dir) const
  { return boundary_side2id_[dir][is_hi]; }
  

  void build_halo(
      int_t num_dims,
      bool with_corners,
      const std::vector<int_t> & block_offsets,
      const std::vector<long_t> & cell_dist,
      std::vector<entity_info_t> & shared_info,
      std::vector<entity_info_t> & ghost_info);

  void pack(int_t id, std::vector<byte_t>& buf, bool as_ghost) const override;
  void unpack(const byte_t * & buf, bool as_ghost) override;
  
  bool is_structured() const override { return true; }
  
  const auto & dims() const { return dims_; }
  auto dim(int_t i) const { return dims_[i]; }
  
  void build_connectivity() override;
  void build_connectivity(int_t a, int_t b) override;
  
  void build_neighbors() override;
  void build_neighbors(int_t dim, int_t thru) override;
  
  void _build_geometry(bool with_face_geom, bool with_cell_geom) override;

  void request_face_geometry() override {}
  void request_cell_geometry() override {}

  using mesh_block_t::output;
#ifdef HAVE_EXODUS
  void output(
      exo_writer_t & exo,
      const std::string & lbl,
      const std::vector<side_set_info_t> & side_sets) const override;
#endif

  void output(vtk_writer_t & vtk) const override;
  
  bool find_cell(const real_t * x, int_t & cell_id) override;
  
  int_t cell_region_id(int_t) const override
  { return region_; }

  void add_intrablock_face(
    const face_iterator_t & f,
    const int_t dir,
    const int_t owner_id,
    const int_t ghost_id);

  int_t add_ghost_entity(long_t gid, int_t dim);
  int_t add_ghost_entity(const int_t * begin, const int_t * end, int_t dim);

  std::vector<real_t> coordinates() const;

  std::vector<std::unique_ptr<mesh_struct_block_t>>
  refine(
      const amr_flags_t & flags,
      const range_ijk_t * ranges,
      int_t nranges ) const;
  
  static std::unique_ptr<mesh_struct_block_t> coarsen(
      int_t split,
      mesh_struct_block_t const * const * blocks);
  
  bool has_neighbor_at(sector_t sec) const;
  const std::vector<struct_neighbor_t> & neighbors_at(sector_t sec) const;
  
  const sector_map_t<std::vector<struct_neighbor_t>> & neighbors() const;
  const std::vector<struct_shared_t> & shared() const;

  void link_intrablock_cells();
};

} // namespace

#endif
