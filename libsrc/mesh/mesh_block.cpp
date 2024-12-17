#include "mesh_block.hpp"
#include "mesh_node.hpp"
#include "mesh_boundary.hpp"
#include "mesh_struct_block.hpp"
#include "mesh_unstruct_block.hpp"

#include "comm/comm_queue.hpp"
#include "io/csv_writer.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Cast as different types
////////////////////////////////////////////////////////////////////////////////
const mesh_struct_block_t * mesh_block_t::as_struct() const
{ return dynamic_cast<const mesh_struct_block_t*>(this); }

mesh_struct_block_t * mesh_block_t::as_struct()
{ return dynamic_cast<mesh_struct_block_t*>(this); }

const mesh_unstruct_block_t * mesh_block_t::as_unstruct() const
{ return dynamic_cast<const mesh_unstruct_block_t*>(this); }

mesh_unstruct_block_t * mesh_block_t::as_unstruct()
{ return dynamic_cast<mesh_unstruct_block_t*>(this); }
  

////////////////////////////////////////////////////////////////////////////////
/// Get id
////////////////////////////////////////////////////////////////////////////////
int_t mesh_block_t::id() const
{ return node_ ? node_->id() : -1; }
  
////////////////////////////////////////////////////////////////////////////////
/// resize block data
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::build_bounding_box() 
{ 
  for (int_t i=0; i<num_dims_; ++i) {
    bounding_box_.lo[i] = consts::real_max;
    bounding_box_.hi[i] = consts::real_min;
  }

  auto num_verts = num_owned_vertices();
  for (int_t i=0, pos=0; i<num_verts; ++i) {
    for (int_t d=0; d<num_dims_; ++d, ++pos) {
      bounding_box_.lo[d] = std::min(vertex_coords_[pos], bounding_box_.lo[d]);
      bounding_box_.hi[d] = std::max(vertex_coords_[pos], bounding_box_.hi[d]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// queue geometry exchange
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::queue_geometry_exchange(comm_queue_block_t & q)
{
  if (cell_centroids_.size())
    q.add(cell_centroids_);
  if (cell_volumes_.size())
    q.add(cell_volumes_);
}

////////////////////////////////////////////////////////////////////////////////
/// request neighbor info
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::build_neighbors(const std::vector<std::pair<int_t, int_t>> & neigh)
{
  for (auto n : neigh) {

    auto dim = n.first;
    auto through = n.second;
    
    if (dim > num_dims_) {
      THROW_ERROR(
          "Neighbor info for 0 <= dim <= " << num_dims_ 
          << ", but you are asking for " << dim);
    }
    if (through > num_dims_) {
      THROW_ERROR(
          "Neighbor info through 0 <= dim <= " << num_dims_ 
          << ", but you are asking for neighbor info through " << dim);
    }
    if (through == dim) {
      THROW_ERROR(
          "Neighbor info through 0 <= dim <= " << num_dims_ 
          << " cannot be found throuogh " << dim);
    }

    desired_neighbors_.emplace(dim, through);
  }
  
  _build_neighbors(neigh);
  
}

////////////////////////////////////////////////////////////////////////////////
/// request connecitivy
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::build_connectivity(const std::vector<std::pair<int_t, int_t>> & conn)
{
  for (auto c : conn) {

    auto a = c.first;
    auto b = c.second;
  
    if (a < 0 || a > num_dims_ || b < 0 || b > num_dims_) {
      THROW_ERROR(
          "Connectivity indices must be 0 <= i <= " << num_dims_ 
          << ", but you are asking for " << a << " -> " << b);
    }
    desired_connectivity_.emplace(a, b); 
  }

  _build_connectivity(conn);
}

////////////////////////////////////////////////////////////////////////////////
/// build connecitivy
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::prune_connectivity()
{
  for (auto it=connectivity_.begin(); it!=connectivity_.end(); ) {
    if (desired_connectivity_.count(it->first)) {
      ++it;
    }
    else {
      it = connectivity_.erase(it);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// get connecitivy
////////////////////////////////////////////////////////////////////////////////
crs_t<int_t,int_t> & mesh_block_t::connectivity(int_t a, int_t b) 
{ 
  auto it = connectivity_.find({a, b});
  if (it == connectivity_.end()) {
    THROW_ERROR("Asking for connectivity " << a << " -> " << b
        << " but this was not generated.  Make sure you request"
        << " this pair beforing calling mesh.initiaize()");
  }
  return it->second;
}

const crs_t<int_t,int_t> & mesh_block_t::connectivity(int_t a, int_t b) const 
{ 
  auto it = connectivity_.find({a, b});
  if (it == connectivity_.end()) {
    THROW_ERROR("Asking for connectivity " << a << " -> " << b
        << " but this was not generated.  Make sure you request"
        << " this pair beforing calling mesh.initiaize()");
  }
  return it->second;
}

////////////////////////////////////////////////////////////////////////////////
/// get neighbors
////////////////////////////////////////////////////////////////////////////////
crs_t<int_t,int_t> & mesh_block_t::neighbors(int_t a, int_t b) 
{ 
  auto it = neighbors_.find({a, b});
  if (it == neighbors_.end()) {
    THROW_ERROR("Asking for neighbors " << a << " -> " << b
        << " but this was not generated.  Make sure you request"
        << " this pair beforing calling mesh.initiaize()");
  }
  return it->second;
}

const crs_t<int_t,int_t> & mesh_block_t::neighbors(int_t a, int_t b) const 
{ 
  auto it = neighbors_.find({a, b});
  if (it == neighbors_.end()) {
    THROW_ERROR("Asking for neighbors " << a << " -> " << b
        << " but this was not generated.  Make sure you request"
        << " this pair beforing calling mesh.initiaize()");
  }
  return it->second;
}


////////////////////////////////////////////////////////////////////////////////
/// build geometry
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::build_geometry(bool with_face_geom, bool with_cell_geom)
{

  if (with_face_geom) {
    face_areas_.clear();
    face_centroids_.clear();
    face_normals_.clear();

    auto nface = num_all_faces();
    face_areas_.resize(nface);
    face_centroids_.resize(nface * num_dims_);
    face_normals_.resize(nface * num_dims_);
  }
  
  if (with_cell_geom) {
    cell_volumes_.clear();
    cell_centroids_.clear();

    auto ncell = num_all_cells();
    cell_volumes_.resize(ncell);
    cell_centroids_.resize(ncell * num_dims_);
  }

  _build_geometry(with_face_geom, with_cell_geom);
}


////////////////////////////////////////////////////////////////////////////////
/// csv outputt
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::output(csv_writer_t & csv)  const
{
  const char * dim_names[] = {"x", "y", "z"};

  csv.comment("Vertex coordinates");;

  for (int_t i=0; i<num_dims_; ++i)
    csv << dim_names[i];

	csv.new_row();

  auto nvert = num_owned_vertices();
  for (int_t p=0; p<nvert; ++p) {
    for (int_t i=0; i<num_dims_; ++i)
      csv << vertex_coords_[p*num_dims_ + i];
    csv.new_row();
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Some general finishing touches
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::map_global_to_local_ids()
{
  auto n = global_ids_.size();
  global_id_map_.clear();
  global_id_map_.resize(n);

  for (size_t i=0; i<n; ++i) {
    const auto & ids = global_ids_[i];
    auto & id_map = global_id_map_[i];
    id_map.clear();
    auto nents = std::min<size_t>(num_owned_entities(i), ids.size());
    for (size_t i=0; i<nents; ++i)
      id_map.emplace( ids[i], i);
  }
}

////////////////////////////////////////////////////////////////////////////////
/// find a cell by its gid
////////////////////////////////////////////////////////////////////////////////
bool mesh_block_t::find_cell(long_t gid, int_t & cell_id)
{
  const auto & id_map = global_id_map_[num_dims_];
  auto it = id_map.find(gid);
  if (it != id_map.end()) {
    cell_id = it->second;
    return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Compute aa global numbering for a disjoint idex space
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::number_disjoint(int_t dim, long_t start)
{
  auto & gids = global_ids(dim);
  gids.clear();
  gids.resize(num_owned_entities(dim));
  std::iota(gids.begin(), gids.end(), start);
}
  
////////////////////////////////////////////////////////////////////////////////
/// Extract boundaries
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::extract_boundaries(
    std::vector<std::unique_ptr<mesh_boundary_t>> & boundaries)
{
  int_t num_bnd = boundaries.size();

  std::vector<boundary_face_t> faces;
  std::vector<boundary_face_t> filtered_faces;

  for (int_t i=0; i<num_bnd; ++i) {
    auto & bnd = boundaries[i];

    auto has_where = bnd->has_where();
      
    faces.clear();
    extract_boundary(bnd->id(), faces);

    const auto & face_vertices = face_centroids_;

    if (has_where && face_vertices.empty())
      THROW_ERROR(
        "You need to build face geometry for coordinate-based " <<
        "boundary filtering"
      );

    // has filter function
    if (has_where) {
      // get function
      auto where_fun = bnd->where();
      // reserve storage
      filtered_faces.clear();
      filtered_faces.reserve(faces.size());
      std::vector<real_t> coords(num_dims_);
      // ffilter out facec ids
      for (auto f : faces) {
        auto start = &face_vertices[f.id*num_dims_];
        coords.assign(start, start + num_dims_);
        if (where_fun(coords))
          filtered_faces.emplace_back(f);
      }
      // set results
      if (!filtered_faces.empty())
        boundary_faces_[i] = filtered_faces;
    }

    // entire boundary 
    else {
      if (!faces.empty())
        boundary_faces_[i] = faces;
    }

  } // boundaries

}

////////////////////////////////////////////////////////////////////////////////
/// Copy metadata
////////////////////////////////////////////////////////////////////////////////
void mesh_block_t::copy_metadata(const mesh_block_t & block) 
{
  desired_neighbors_ = block.desired_neighbors_;
  desired_connectivity_ = block.desired_connectivity_;
}

} // naemspace prl
