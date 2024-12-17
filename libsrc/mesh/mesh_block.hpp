#ifndef MESH_BLOCK_HPP
#define MESH_BLOCK_HPP

#include "config.hpp"
#include "face_data.hpp"
#include "io/side_set.hpp"
#include "geom/bounding_box.hpp"
#include "utils/array_ref.hpp"
#include "utils/crs.hpp"
#include "utils/errors.hpp"

#include <map>
#include <memory>
#include <set>

namespace prl {

class mesh_boundary_t;
class mesh_node_t;
class mesh_struct_block_t;
class mesh_unstruct_block_t;
struct comm_map_block_t;
class comm_queue_block_t;
class csv_writer_t;
class exo_writer_t;
class vtk_writer_t;

////////////////////////////////////////////////////////////////////////////////
/// Base mesh block class
////////////////////////////////////////////////////////////////////////////////
class mesh_block_t {
protected:

  int_t num_dims_ = 0;
  mesh_node_t * node_ = nullptr;

public:
  
  bounding_box_t bounding_box_;
  
  std::map< int_t, std::vector<boundary_face_t> > boundary_faces_;
  
  std::vector<real_t> vertex_coords_;
  
  std::vector<real_t> face_centroids_;
  std::vector<real_t> face_normals_;
  std::vector<real_t> face_areas_;

  std::vector<real_t> cell_centroids_;
  std::vector<real_t> cell_volumes_;

  std::vector<int_t> num_all_entities_;
  std::vector<int_t> num_owned_entities_;
  
  std::vector<std::vector<long_t>> global_ids_;
  std::vector<std::map<long_t, int_t>> global_id_map_;

  std::set<std::pair<int_t, int_t>> desired_connectivity_;
  std::map< std::pair<int_t, int_t>, crs_t<int_t, int_t> > connectivity_;
  
  std::set<std::pair<int_t, int_t>> desired_neighbors_;
  std::map< std::pair<int_t, int_t>, crs_t<int_t, int_t> > neighbors_;

  mesh_block_t(int_t num_dims) :
    num_dims_(num_dims)
  { }
  
  const mesh_struct_block_t * as_struct() const;
  mesh_struct_block_t * as_struct();
  
  const mesh_unstruct_block_t * as_unstruct() const;
  mesh_unstruct_block_t * as_unstruct();
  
  int_t id() const;

  auto node() { return node_; }
  void set_node(mesh_node_t * n) { node_ = n; }

  auto num_dims() const { return num_dims_; }

  // dim = 0 :vertex. 
  const std::vector<long_t> & global_ids(size_t dim) const
  { return global_ids_[dim]; }

  std::vector<long_t> & global_ids(size_t dim)
  { 
    if (dim >= global_ids_.size()) global_ids_.resize(dim+1);
    return global_ids_[dim];
  }

  bool has_boundary(int_t i) { return boundary_faces_.count(i); }

  void extract_boundaries(std::vector<std::unique_ptr<mesh_boundary_t>> & boundaries);

  virtual void extract_boundary(
    int_t boundary_id,
    std::vector<boundary_face_t> & faces) = 0;

  void build_neighbors(const std::vector<std::pair<int_t,int_t>> & conn);
  virtual void _build_neighbors(const std::vector<std::pair<int_t,int_t>> & conn) = 0;
  virtual void _build_neighbors(int_t a, int_t  b) = 0;
  
  virtual void pack(int_t id, std::vector<byte_t>& buf, bool as_ghost) const = 0;
  virtual void unpack(const byte_t * & buf, bool as_ghost) = 0;
  
  void map_global_to_local_ids();
  void build_bounding_box();

  void queue_geometry_exchange(comm_queue_block_t & q);

  void prune_connectivity();

  void build_connectivity(const std::vector<std::pair<int_t,int_t>> & conn);
  virtual void _build_connectivity(const std::vector<std::pair<int_t,int_t>> & conn) = 0;
  virtual void _build_connectivity(int_t a, int_t  b) = 0;
  virtual void _build_connectivity_post() {}


  crs_t<int_t, int_t> & connectivity(int_t a, int_t b);
  const crs_t<int_t, int_t> & connectivity(int_t a, int_t b) const;
  
  crs_t<int_t, int_t> & neighbors(int_t from, int_t thru);
  const crs_t<int_t, int_t> & neighbors(int_t from, int_t thru) const;

  void build_geometry(bool with_face_geom, bool with_cell_geom);
  virtual void _build_geometry(bool with_face_geom, bool with_cell_geom) = 0;
  
  int_t num_all_entities(size_t dim) const
  { return num_all_entities_.at(dim); }
  
  int_t & num_all_entities(size_t dim)
  {
    if (dim >= num_all_entities_.size()) num_all_entities_.resize(dim+1);
    return num_all_entities_.at(dim);
  }

  int_t num_owned_entities(size_t dim) const
  { return num_owned_entities_.at(dim); }
  
  int_t & num_owned_entities(size_t dim)
  {
    if (dim >= num_owned_entities_.size()) num_owned_entities_.resize(dim+1);
    return num_owned_entities_.at(dim);
  }

  int_t num_owned_cells() const
  { return num_owned_entities(num_dims_); }
  int_t num_all_cells() const
  { return num_all_entities(num_dims_); }
  
  int_t num_all_faces() const
  { return num_all_entities(num_dims_-1); }
  int_t num_owned_faces() const
  { return num_owned_entities(num_dims_-1); }
  
  int_t num_all_vertices() const
  { return num_all_entities(0); }
  int_t num_owned_vertices() const
  { return num_owned_entities(0); }

  virtual void output(csv_writer_t & csv) const;
#ifdef HAVE_EXODUS
  virtual void output(
      exo_writer_t & exo,
      const std::string & label,
      const std::vector<side_set_info_t> & side_sets = {}) const = 0;
#endif
  virtual void output(vtk_writer_t & vtk) const = 0;
  
  virtual bool is_structured() const = 0;

  virtual bool find_cell(const real_t * x, int_t & cell_id) = 0;
  bool find_cell(long_t gid, int_t & cell_id);

  virtual int_t cell_region_id(int_t id) const = 0;

  void number_disjoint(int_t dim, long_t start);

  void copy_metadata(const mesh_block_t & block);
  
  virtual ~mesh_block_t() {}
};

} // namespace

#endif
