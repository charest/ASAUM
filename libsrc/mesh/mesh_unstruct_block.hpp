#ifndef UNSTRUCTURED_BLOCK_HPP
#define UNSTRUCTURED_BLOCK_HPP

#include "mesh_block.hpp"
#include "io/shape.hpp"
#include "utils/crs.hpp"

#include <map>
#include <memory>
#include <vector>

namespace prl {

class mesh_boundary_t;

////////////////////////////////////////////////////////////////////////////////
/// Structured mesh block
////////////////////////////////////////////////////////////////////////////////
class mesh_unstruct_block_t : public mesh_block_t {
    
public:
  
  std::map<int_t, std::vector<long_t>> local2global_;
  std::map<int_t, std::map<long_t, int_t>> global2local_;
 
  std::vector<shape_t> cell_type_;
  std::vector<int_t> cell_region_;

  std::vector<int_t> face_owner_; // 3d only
  std::map<std::vector<int_t>, int_t> sorted_face_vertices_; // 3d only

  std::vector<int_t> regions_;
  std::vector<int_t> region_cell_offsets_;

  std::vector<int_t> side_set_offsets_;
  std::vector<int_t> side_sets_;

  std::vector<int_t> side_id_;
  std::vector<int_t> side_cell_;
  std::vector<int_t> side_cell_face_;

  std::map<long_t, std::vector<int_t>> cell_sides_;
  
  std::map<int_t, std::vector<int_t>> boundary_side2id_;
  std::map<int_t, int_t> boundary_id2side_;

  mesh_unstruct_block_t(int_t num_dims);

  void cell_midpoints(std::vector<real_t> &);

  void pack(int_t id, std::vector<byte_t>& buf, bool as_ghost) const override;
  void unpack(const byte_t * & buf, bool as_ghost) override;

  void erase(const std::vector<size_t> & ids);
  
  void add_region(int_t id, int_t num);
  
  void build_cell_sides();

  void group_regions();
  void group_sides();

  void sort_face_vertices();
  
  void reorder_cells(const std::vector<int_t> & order);
  void reorder_sides(const std::vector<int_t> & order);

  bool is_structured() const override { return false; }
  
  void transpose(int_t from, int_t to);
  void build(int_t dim);
  void intersect(int_t from, int_t to, int_t through);
  
  void build_connectivity(int_t a, int_t b) override;
  void build_connectivity() override;
  
  void build_neighbors() override;
  void build_neighbors(int_t dim, int_t thru) override;
  
  void post_connectivity() override;
  
  void _build_geometry(bool with_face_geom, bool with_cell_geom) override;
  
  void request_face_geometry() override;
  void request_cell_geometry() override;
  
  using mesh_block_t::output;

#ifdef HAVE_EXODUS
  void output(
      exo_writer_t & exo,
      const std::string & lbl,
      const std::vector<side_set_info_t> & side_sets) const override;
#endif
  void output(vtk_writer_t & vtk) const override;

  bool find_cell(const real_t * x, int_t & cell_id) override;
  
  int_t cell_region_id(int_t id) const override
  { return cell_region_[id]; }

  void set_num_vertices();
  void set_num_cells();
  void set_num_entities(int_t dim);
  
  void extract_boundary(int_t boundary_id, std::vector<boundary_face_t> & faces) override;

  void add_boundary_mapping(int_t side_id, int_t bid);
  void add_boundary_mappings(int_t side_id, const std::vector<int_t> labels);

  bool find_boundary_from_id(int_t bid, int_t & side_id) const;

  void install_boundaries(
      const std::vector<std::unique_ptr<mesh_boundary_t>> & boundaries);
};

} // namespace

#endif
