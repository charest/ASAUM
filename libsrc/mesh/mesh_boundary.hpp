#ifndef MESH_BOUNDARY_HPP
#define MESH_BOUNDARY_HPP

#include "common/boundary.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Mesh Boundary info
////////////////////////////////////////////////////////////////////////////////
class mesh_boundary_t : public boundary_t {
public:
  mesh_boundary_t(
      int_t id,
      const std::string & name)
    : boundary_t(id, name)
  {}
  virtual ~mesh_boundary_t() = default;
};

////////////////////////////////////////////////////////////////////////////////
/// Mesh Boundary info
////////////////////////////////////////////////////////////////////////////////
class mesh_struct_boundary_t : public mesh_boundary_t {
  int_t side_id_ = -1;
  std::vector<std::vector<int_t>> patches_;
public:
  mesh_struct_boundary_t(
      int_t id,
      const std::string & name,
      const std::vector<std::vector<int_t>> & vs) 
    : mesh_boundary_t(id, name), patches_(vs)
  {}
  const auto & patches() const { return patches_; }
};

////////////////////////////////////////////////////////////////////////////////
/// Mesh Boundary info
////////////////////////////////////////////////////////////////////////////////
class mesh_unstruct_boundary_t : public mesh_boundary_t {
  int_t side_id_ = -1;
public:
  mesh_unstruct_boundary_t(
      int_t id,
      const std::string & name,
      int_t side_id)
    : mesh_boundary_t(id, name), side_id_(side_id)
  {}
  const auto side_id() const { return side_id_; }
};


////////////////////////////////////////////////////////////////////////////////
/// Boundary helper
////////////////////////////////////////////////////////////////////////////////
std::unique_ptr<mesh_boundary_t> make_mesh_boundary(
    lua_result_t input,
    mpi_comm_t & comm,
    int_t pos,
    bool is_structured);

} // namespace

#endif
