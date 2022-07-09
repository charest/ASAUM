#include "mesh_boundary.hpp"

#include "comm/mpi_comm.hpp"
#include "lua/lua_utils.hpp"

#include <algorithm>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Boundary helper
////////////////////////////////////////////////////////////////////////////////
std::unique_ptr<mesh_boundary_t> make_mesh_boundary(
    lua_result_t input,
    mpi_comm_t & comm,
    int_t pos,
    bool is_structured)
{

  auto is_root = comm.is_root();

  //--- Name
  auto res = as_scalar<std::string>(input, "name", is_root);
  if (!res.exists) comm.exit(consts::failure);
  const auto & name = res.value;

  std::unique_ptr<mesh_boundary_t> bc;


  //--- Structured reads face connectivity
  if (is_structured) {

    if (!validate(input, "connect", is_root))
      comm.exit(consts::failure);

    auto conn_input = input["connect"];
    auto nconn = conn_input.size();
    if (nconn < 1) {
      if (is_root) {
        std::cout << "Need at least one list of vertices for boundary face";
        std::cout << " connectivity." << std::endl;
      }
      comm.exit(consts::failure);
    }
  
    std::vector< std::vector<int_t> > connectivity;

    for (decltype(nconn) i=0; i<nconn; ++i) {
      auto conn = conn_input[i+1].as<std::vector<int_t>>();
      if (!conn.size()) {
        if (is_root) {
          std::cout << "Boundary connectivity entry must have at least one value.";
          std::cout << std::endl;
        }
        comm.exit(consts::failure);
      }
      std::sort(conn.begin(), conn.end());
      connectivity.emplace_back(conn);
    }

    bc = std::make_unique<mesh_struct_boundary_t>(pos, name, connectivity);
  }

  //--- Unstructured reads ids
  else {
    auto res = as_scalar<int_t>(input, "id", is_root); 
    if (!res.exists)
      comm.exit(consts::failure);
    bc = std::make_unique<mesh_unstruct_boundary_t>(pos, name, res.value);
  }
  

  //--- Apply filtering if requested
  auto where_input = input["where"];
  if (!where_input.empty()) {
    auto where_fun = [=](const std::vector<real_t> & x) -> bool 
    {
      auto ret = where_input(x);
      if (!ret.is_representable_as<bool>())
        THROW_ERROR("Could not represent the result of 'where' as a boolean!");
      return ret.as<bool>();
    };
    bc->where( where_fun );
  }

  return bc;
}

} // namespace
