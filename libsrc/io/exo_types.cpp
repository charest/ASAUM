#include "exo_types.hpp"

namespace prl {

/* triangle */
constexpr int tri_sides = 3;
constexpr int tri_size = 3;
static int tri_table[tri_sides][tri_size] = {
  /*   1        2        3                                             side   */
  {1,2,4}, {2,3,5}, {3,1,6}                                       /* nodes  */
};

/* quad */
constexpr int quad_sides = 4;
constexpr int quad_size = 3;
static int quad_table[quad_sides][quad_size] = {
  /*   1        2        3        4                                    side   */
  {1,2,5}, {2,3,6}, {3,4,7}, {4,1,8}                              /* nodes  */
};

/* tetra */
constexpr int tetra_sides = 4;
constexpr int tetra_size = 6;
static int tetra_table[tetra_sides][tetra_size] = {
  /*      1              2               3               4            side   */
  {1,2,4,5,9,8}, {2,3,4,6,10,9}, {1,4,3,8,10,7}, {1,3,2,7,6,5}   /* nodes  */
};

/* hex */
constexpr int hex_sides = 6;
constexpr int hex_size = 9;
static int hex_table[hex_sides][hex_size] = {
  /*         1                        2                                side   */
  {1,2,6,5,9,14,17,13,26}, {2,3,7,6,10,15,18,14,25},              /* nodes  */
  /*         3                        4                                side   */
  {3,4,8,7,11,16,19,15,27}, {1,5,8,4,13,20,16,12,24},             /* nodes  */
  /*         5                        6                                side   */
  {1,4,3,2,12,11,10,9,22},  {5,6,7,8,17,18,19,20,23}              /* nodes  */
};


/* shell */
constexpr int shell_sides = 6;
constexpr int shell_size = 8;
static int shell_table[shell_sides][shell_size] = {
  /*        1                  2                                       side   */
  {1,2,3,4,5,6,7,8}, {1,4,3,2,8,7,6,5} ,                          /* nodes  */
  /*        3                  4                                       side   */
  {1,2,5,0,0,0,0,0}, {2,3,6,0,0,0,0,0} ,                          /* nodes  */
  /*        5                  6                                       side   */
  {3,4,7,0,0,0,0,0}, {4,1,8,0,0,0,0,0}                            /* nodes  */
};

void elem_side_table(
    int elem_type,
    int *& side_vert_map,
    int & table_len,
    int & num_sides)
{
  side_vert_map = nullptr;
  table_len = 0;
  num_sides = 0;

  switch (elem_type) {
  case (EX_EL_TRIANGLE):
    table_len = tri_size;
    num_sides = tri_sides;
    side_vert_map = &tri_table[0][0];
    return;
  case (EX_EL_QUAD):
    table_len = quad_size;
    num_sides = quad_sides;
    side_vert_map = &quad_table[0][0];
    return;
  case (EX_EL_SHELL):
    table_len = shell_size;
    num_sides = shell_sides;
    side_vert_map = &shell_table[0][0];
    return;
  case (EX_EL_TETRA):
    table_len = tetra_size;
    num_sides = tetra_sides;
    side_vert_map = &tetra_table[0][0];
    return;
  case (EX_EL_HEX):
    table_len  = hex_size;
    num_sides = hex_sides;
    side_vert_map = &hex_table[0][0];
    return;
  default:
    return;
  };
}
  
} // namespace
