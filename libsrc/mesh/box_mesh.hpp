#ifndef BOX_MESH_HPP
#define BOX_MESH_HPP

#include "mesh.hpp"
#include "ijk/mesh_ijk.hpp"
#include "utils/result.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Box mesh class
////////////////////////////////////////////////////////////////////////////////
class box_mesh_t : public mesh_t {
  
  std::unique_ptr<mesh_ijk_t> initial_mesh_;
  
public:

  box_mesh_t(
      lua_result_t input,
      mpi_comm_t & comm,
      const std::vector<int_t> & parts);
  
  bool is_structured() const override { return true; }

  void load() override;

  void number() override;
  void build_halo(std::vector<comm_map_block_t> & comm_maps) override;

  void number_vertices();

};

}

#endif
