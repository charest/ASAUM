#ifndef DRIVER_HPP
#define DRIVER_HPP

#include "comm/comm_map.hpp" 

#include "io/csv_writer.hpp"
#ifdef HAVE_EXODUS
#include "io/exo_writer.hpp"
#endif
#include "io/vtk_writer.hpp"
#include "lua/lua_utils.hpp"
#include "mesh/mesh.hpp"

#include <memory>
#include <string>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Unstructured driver
////////////////////////////////////////////////////////////////////////////////
struct driver_t {

  lua_t input_;

  std::vector<int_t> refine_;
  
  /// counters
  int_t output_counter_ = 0;
  real_t clock_start_ = 0;
  real_t solution_time_ = 0;

  /// output
  std::string output_prefix_;

  /// amr
  int_t amr_levels_ = 5;
  int_t initial_amr_ = 0;

  // packages
  std::unique_ptr<mesh_t> mesh_;

  // commmunication data
  mpi_comm_t & comm_;

  // io data
  std::unique_ptr<csv_writer_t> csv_;
#ifdef HAVE_EXODUS
  std::unique_ptr<exo_writer_t> exo_;
#endif
  std::unique_ptr<vtk_writer_t> vtk_;

  /// constructor
  driver_t(
      lua_t input,
      mpi_comm_t & comm);
  
  void initialize();
  
  void amr();
  void periodic_amr();
  void initial_amr();

  void refine(
    const amr_flags_t & flags,
    mesh_block_t * parent,
    mesh_block_t ** children,
    int_t nchild);
  void coarsen(
    const amr_flags_t & flags,
    mesh_block_t * parent,
    mesh_block_t ** children,
    int_t nchild);

  void output();
  
  void create_directory(const std::string & dirname);

  void run();

};

} // namespace

#endif
