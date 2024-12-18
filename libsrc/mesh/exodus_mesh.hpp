#ifndef EXODUS_MESH_HPP
#define EXODUS_MESH_HPP

#include "mesh.hpp"
#include "partition.hpp"
#include "io/exo_types.hpp"
#include "io/side_set.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Exodus mesh class
////////////////////////////////////////////////////////////////////////////////
class exodus_mesh_t : public mesh_t {
  

  std::string filename_;

  distribution_t distribution_ = distribution_t::sequential;
  part_alg_t part_alg_ = part_alg_t::parmetis_kway;
  bool repart_ = false;

public:
  
  using side_list = std::vector<side_set_info_t>;
  side_list side_sets_;
  std::map<int_t, side_list::iterator> side_id_map_;
  std::map<std::string, side_list::iterator> side_name_map_;

  exodus_mesh_t(
      lua_result_t input,
      mpi_comm_t & comm,
      const std::vector<int_t> & parts,
      const std::string & filename,
      distribution_t distribution,
      part_alg_t part_alg,
      bool repart);
  
  bool is_structured() const override { return false; }
  
  void read_and_partition();

  void read(
    mpi_comm_t & comm,
    bool serial,
    const std::vector<std::string> & filenames,
    const std::vector<int_t> & part_ids);

  void add_side_set(const side_set_info_t &);
  
  void post_migration(mpi_comm_t & comm);

  void number() override;
  void _build_halo(
      int_t num_ghost,
      bool with_corners,
      std::vector<comm_map_block_t> & comm_maps) override;

  void output(exo_writer_t & exo, int_t blk, const std::string & lbl) const override
  { blocks_[blk]->output(exo, lbl, side_sets_); }

  void migrate(
      const std::vector<int_t> & cell_parts,
      const std::vector<int_t> & part_dist,
      const std::vector<int_t> & part_ids);

  void graph(
      int_t min_connections,
      std::vector<size_t> & distribution,
      crs_t<size_t, int_t> & graph);
};

} // namespace

#endif
