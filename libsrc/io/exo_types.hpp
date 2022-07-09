#ifndef EXO_BLOCK_HPP
#define EXO_BLOCK_HPP

#include "config.hpp"

#include <exodusII.h>

#include <map>
#include <string>
#include <vector>

namespace prl {

template<typename U>
struct exo_block_stats_t {
  U num_elem_this_blk = 0;
  U num_faces_per_elem = 0;
  U num_edges_per_elem = 0;
  U num_nodes_per_elem = 0;
  U num_attr = 0;
  char elem_type[MAX_STR_LENGTH];
};

enum ex_element_type {
  EX_EL_UNK         =  -1,     /**< unknown entity */
  EX_EL_NULL_ELEMENT=   0,
  EX_EL_TRIANGLE    =   1,     /**< Triangle entity */
  EX_EL_QUAD        =   2,     /**< Quad entity */
  EX_EL_HEX         =   3,     /**< Hex entity */
  EX_EL_WEDGE       =   4,     /**< Wedge entity */
  EX_EL_TETRA       =   5,     /**< Tetra entity */
  EX_EL_TRUSS       =   6,     /**< Truss entity */
  EX_EL_BEAM        =   7,     /**< Beam entity */
  EX_EL_SHELL       =   8,     /**< Shell entity */
  EX_EL_SPHERE      =   9,     /**< Sphere entity */
  EX_EL_CIRCLE      =  10,     /**< Circle entity */
  EX_EL_TRISHELL    =  11,     /**< Triangular Shell entity */
  EX_EL_PYRAMID     =  12      /**< Pyramid entity */
};

struct ex_elem_blk_parm
{
  char elem_type[33];
  int64_t elem_blk_id;
  int64_t num_elem_in_blk;
  int num_nodes_per_elem;
  int num_sides;
  int num_nodes_per_side[6];
  int num_attr;
  int64_t elem_ctr;
  ex_element_type elem_type_val;
};
  
//! side set stats
template<typename U>
struct exo_side_set_stats_t {
  U num_side_in_set;
  U num_dist_fact_in_set;
  U side_set_node_list_len;
  U tot_num_ss_elem;
};
  

void elem_side_table(
    int elem_type,
    int *& side_vert_map,
    int & table_len,
    int & num_sides);

} // namespace

#endif
