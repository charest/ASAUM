#ifndef PARTITION_HPP
#define PARTITION_HPP

#include "config.hpp"

#include <stddef.h>
#include <vector>

namespace prl {

class mpi_comm_t;

template<typename T, typename U>
struct crs_t;

// the different partiton algorithms
enum class part_alg_t {
  parmetis_kway,
  parmetis_geom,
  parmetis_geomkway,
  metis,
  naive
};

// make a graph
void graph(
    mpi_comm_t & comm,
    const std::vector<int_t> block_offsets,
    const std::vector< const crs_t<int_t, int_t> * > & cell2vert,
    const std::vector< const std::vector<long_t> * > & vert_local2global,
    int_t min_connections,
    std::vector<size_t> & distribution,
    crs_t<size_t, int_t> & graph);

/// partition a mesh
void partition(
    mpi_comm_t & comm,
    const std::vector<int_t> & block_offsets,
    const std::vector<size_t> & distribution,
    const crs_t<size_t, int_t> & graph,
    part_alg_t part_alg,
    const std::vector<real_t> & cell_xyz,
    int_t tot_parts,
    std::vector<int_t> & parts);

} // namespace

#endif
