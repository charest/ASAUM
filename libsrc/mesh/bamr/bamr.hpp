#ifndef AMR_HPP
#define AMR_HPP

#include "config.hpp"
#include "bamr_block.hpp"
#include "comm/mpi_comm.hpp"
#include "mesh/mesh.hpp"

#include <vector>

namespace prl {

struct bamr_flags_t;

////////////////////////////////////////////////////////////////////////////////
// Main AMR class
////////////////////////////////////////////////////////////////////////////////
class bamr_t {
public:

  using refine_function_t = std::function<void(int_t, bamr_block_t**, int_t)>;

private:

  int_t max_level_ = 10;

  int_t num_dims_ = -1;

  mpi_comm_t & comm_;
  
  std::vector<bamr_block_t> blocks_;

public:

  bamr_t(int_t ndims, mpi_comm_t & comm);

  void reserve_blocks(size_t nblocks)
  { blocks_.reserve(nblocks); }

  bamr_block_t & add_block();
  bamr_block_t & add_block(int_t label);

  size_t count_leaves() const;
  
  void get_leaves(std::vector<bamr_block_t*> & list) const;
  std::vector<bamr_block_t*> get_leaves() const;

  void fixup(
      std::vector<bamr_block_t*> & blocks,
      std::vector<bamr_flags_t> & flags);
  
  bool check(
      std::vector<bamr_block_t*> & blocks,
      std::vector<bamr_flags_t> & flags,
      bool first_pass);

  bool inner_check(
      std::vector<bamr_block_t*> & blocks,
      std::vector<bamr_flags_t> & flags,
      bool first_pass);
  
  bool final_check(
      std::vector<bamr_block_t*> & blocks,
      std::vector<bamr_flags_t> & flags);

  void change_connect(
      std::vector<bamr_block_t*> & blocks,
      const std::vector<bamr_flags_t> & flags);

  void refine_tree(
      std::vector<bamr_block_t*> & blocks,
      const std::vector<bamr_flags_t> & flags);
  
  void coarsen_tree(
      std::vector<bamr_flags_t> & flags,
      int_t ndim);

  void amr(const std::vector<bamr_flags_t> & flags, refine_function_t func);
  void refine_uniformly(const bool * flags, refine_function_t func);
  void refine_uniformly(int_t * levels, refine_function_t func);

  void coarsen_uniformly();

  void number();

};

} // namespace

#endif
