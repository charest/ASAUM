#include <exodusII.h>
#include <mpi.h>

#include <cstring>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

#include "math/decomp.hpp"
#include "math/subdivide.hpp"
#include "utils/arguments.hpp"

using int_t = int;
using ulong_t = std::size_t;
using real_t = double;

using namespace prl;

// Command line arguments
option_t<std::vector<ulong_t>> ArgDims(
    "dimensions",
    option_description_t("Box dimensions."),
    required_t::yes);

option_t<std::vector<int_t>> Parts(
    "partitions",
    option_description_t("Number partitions in each direction."));


option_t<std::vector<real_t>> Lower(
      "lower",
      option_description_t("Box lower coordinates."));


option_t<std::vector<real_t>> Upper(
      "upper",
      option_description_t("Box upper coordinates."));


option_t<bool> Verbose(
    "verbose",
    option_description_t("Print extra debug info."));

option_t<bool> LargeInteger(
    "large-integer",
    option_description_t("Use large integer support."));

option_t<std::string> OutputFile(
    "output-file",
    option_description_t("The The output file.."),
    value_description_t("filename"),
    required_t::yes);


/// side set datastructure
struct side_t {
  std::string name;
  std::vector<ulong_t> elem_list;
  std::vector<ulong_t> side_list;
};
 
////////////////////////////////////////////////////////////////////////////////
/// write cell info 
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Container>
auto write_cells(
  int exo_id,
  ex_entity_id blk_id,
  const Container & cell_vertices,
  ulong_t num_verts_per_cell
) {
  
  std::string desc;
  
  switch(num_verts_per_cell) {
    case 2:
      desc = "bar";
      break;
    case 4:
      desc = "quad4";
      break;
    case 8:
      desc = "hex8";
      break;
    default:
      std::cerr
        << "Unknown element type with " << num_verts_per_cell << " verts"
        << std::endl;
      return -1;
  }
  
  auto num_elem = cell_vertices.size() / num_verts_per_cell;
  auto status = ex_put_block( exo_id, EX_ELEM_BLOCK, blk_id, desc.c_str(),
    num_elem, num_verts_per_cell, 0, 0, 0 );
  if (status) return status;

  status = ex_put_name(exo_id, EX_ELEM_BLOCK, blk_id, "element block");
  if (status) return status;

  if ( std::is_same<typename Container::value_type, T>::value )
    status = ex_put_elem_conn( exo_id, blk_id, cell_vertices.data() );
  else {
    std::vector<T> vec(cell_vertices.begin(), cell_vertices.end());
    status = ex_put_elem_conn( exo_id, blk_id, vec.data() );
  }
  
  return status;
}

////////////////////////////////////////////////////////////////////////////////
/// write global id mapping
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Container>
auto write_mapping(
  int exo_id,
  const Container & cell_mapping,
  const Container & vert_mapping
) {
  int status;
  if ( std::is_same<typename Container::value_type, T>::value ) {
    status = ex_put_elem_num_map( exo_id, cell_mapping.data() );
    if (status) return status;
    status = ex_put_node_num_map( exo_id, vert_mapping.data() );
    if (status) return status;
  }
  else {
    std::vector<T> vec(cell_mapping.begin(), cell_mapping.end());
    status = ex_put_elem_num_map( exo_id, vec.data() );
    if (status) return status;
    vec.assign(vert_mapping.begin(), vert_mapping.end());
    status = ex_put_node_num_map( exo_id, vec.data() );
    if (status) return status;
  }
  return status;
}

////////////////////////////////////////////////////////////////////////////////
/// write side set
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Container>
auto write_side_set(
    int ex_id,
    int side_id,
    const Container & elem_list,
    const Container & side_list)
{
  auto status = ex_put_side_set_param(ex_id, side_id, side_list.size(), 0);
  if (status) return status;

  if ( std::is_same<typename Container::value_type, T>::value ) {
    return ex_put_side_set(ex_id, side_id, elem_list.data(), side_list.data());
  }
  else {
    std::vector<T> elem(elem_list.begin(), elem_list.end());
    std::vector<T> side(side_list.begin(), side_list.end());
    return ex_put_side_set(ex_id, side_id, elem.data(), side.data());
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Main driver
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char* argv[] )
{
  
  bool is_verbose = false;
  bool force_int64 = false;

  MPI_Init(&argc, &argv);

  int comm_rank, comm_size;
  MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
  MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );

  //============================================================================
  // parse command line arguments
  auto ret = parse_command_line_arguments(argc, argv, comm_rank==0);
  if (ret > 0) {
    MPI_Finalize();
    return 0;
  }
  else if (ret < 0) {
    MPI_Finalize();
    return 1;
  }

  bool print_help = false;
  
  std::string name = OutputFile;

  std::vector<ulong_t> global_dims = ArgDims;
  int_t num_dims = global_dims.size();

  std::vector<real_t> global_lower(num_dims, -0.5);
  std::vector<real_t> global_upper(num_dims,  0.5);
  
  // check the number of dimensions
  if (num_dims > 3 || num_dims < 1) {
    if (comm_rank == 0) {
      std::cout << "Only one, two to three dimensions supported!" << std::endl;
    }
    print_help = true;
  }

  //============================================================================
  // get the lower and upper coordinates of the bounding box
  if (Lower) {
    std::vector<real_t> args = Lower;
    int_t n = args.size();
    if ( n != num_dims ) {
      if (comm_rank == 0) {
        std::cout
          << "Number of lower box coordinates must match the number of "
          << "dimensions, " << n << " provided, " << num_dims << " expected."
          << std::endl;
      }
      print_help = true;
    }
    global_lower = args;
  }

  if (Upper) {
    std::vector<real_t> args = Upper;
    int_t n = args.size();
    if ( n != num_dims ) {
      if (comm_rank == 0) {
        std::cout
          << "Number of upper box coordinates must match the number of "
          << "dimensions, " << n << " provided, " << num_dims << " expected."
          << std::endl;
      }
      print_help = true;
    }
    global_upper = args;
  }

  if (print_help) {
    MPI_Finalize();
    return -1;
  }

  // get partitioning
  std::vector<ulong_t> block_sizes(num_dims, 1);

  //------------------------------------
  // Partitioning specified
  if (Parts) {
    if (comm_rank==0) std::cout << "Using provided partitioning" << std::endl;
    std::vector<int_t> args = Parts;
    int_t n = args.size();
    if ( n != num_dims ) {
      if ( comm_rank == 0 ) {
        std::cout
          << "Number of partition must match the number of "
          << "dimensions, " << n << " privided, " << num_dims << " expected."
          << std::endl;
      }
      MPI_Finalize();
      return -1;
    }
    ulong_t total_partitions = 1;
    for ( int_t i=0; i<num_dims; ++i ) {
      block_sizes[i] = args[i];
      total_partitions *= args[i];
    }
    if ( total_partitions != static_cast<ulong_t>(comm_size) ) {
      if ( comm_rank == 0 ) {
        std::cout
          << "Total number of partitions does not match the total number of "
          << " mpi ranks" << std::endl;
      }
      MPI_Finalize();
      return -1;
    }
  }
  //------------------------------------
  // figure out the partitioning
  else {
    if (num_dims == 2) {
      block_sizes = decomp(global_dims, comm_size);
    } else if (num_dims == 3) {
      block_sizes = decomp(global_dims, comm_size);
    } else {
      block_sizes[0] = comm_size;
    }

    if (comm_rank==0) {
      std::cout << "Automatically partitioned block sizes: ";
      for ( auto i : block_sizes ) std::cout << i << " ";
      std::cout << std::endl;
    }

  } // partitioning
  //------------------------------------
  

  if (LargeInteger) {
    force_int64 = true;
    if (comm_rank == 0) {
      std::cout << "Forcing 64 bit integers" << std::endl;
    }
  }

  if (Verbose) {
    is_verbose = true;
    if (comm_rank == 0) {
      std::cout << "Adding extra exodus debug." << std::endl;
    }
    ex_opts(EX_ABORT | EX_VERBOSE);
  }

  if (comm_rank == 0) {
    std::cout << num_dims << "-dimensional mesh selected." << std::endl;
  }

  //============================================================================
  // Figureout the dimensions of each local partition
  std::vector<real_t> delta(num_dims);
  for ( int_t i=0; i<num_dims; ++i ) {
    if ( global_upper[i] < global_lower[i] ) {
      if (comm_rank == 0) {
        std::cout
          << "Bounding box is invalid: for i=" << i << ", " << global_lower[i]
          << " should be less than " << global_upper[i]
          << std::endl;
      }
      MPI_Finalize();
      return -1;
    }
    delta[i] =  (global_upper[i] - global_lower[i]) / global_dims[i];
  }

  std::vector<ulong_t> block_ijk(num_dims, comm_rank);
  auto temp = comm_rank;
  for ( int_t i=0; i<num_dims-1; ++i ) {
    block_ijk[i] = temp % block_sizes[i];
    block_ijk[i+1] = (temp - block_ijk[i]) / block_sizes[i];
    temp -= block_ijk[i];
    temp /= block_sizes[i];
  }

  if (is_verbose) {
    std::cout << std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
    for ( int i=0; i<comm_size; ++i )
    {
      if ( comm_rank == i ) {
        std::cout << "Rank " << comm_rank << " has block index: ";
        for ( int_t j=0; j<num_dims; ++j ) std::cout << block_ijk[j] << " ";
        std::cout << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  
  std::vector< std::vector<ulong_t> > block_displ(num_dims);
  for ( int_t i=0; i<num_dims; ++i )
    subdivide( global_dims[i], block_sizes[i], block_displ[i] );

  ulong_t num_verts = 1;
  ulong_t num_cells = 1;
  ulong_t global_verts = 1;
  ulong_t global_cells = 1;
  ulong_t num_verts_per_cell = 1;

  std::vector<ulong_t> local_dims(num_dims);
  std::vector<real_t> local_lower(num_dims);
  std::vector<real_t> local_upper(num_dims);

  for ( int_t i=0; i<num_dims; ++i ) {
    global_verts *= global_dims[i]+1;
    global_cells *= global_dims[i];
    num_verts_per_cell *= 2;
    auto jblock = block_ijk[i];
    const auto & displ = block_displ[i];
    local_dims[i] = displ[jblock+1] - displ[jblock];
    num_cells *= local_dims[i];
    num_verts *= local_dims[i]+1;
    local_lower[i] = global_lower[i] + (displ[jblock])*delta[i] ;
    local_upper[i] = global_lower[i] + (displ[jblock+1])*delta[i];
  }


  MPI_Datatype mpi_type_t;
  if ( sizeof(ulong_t) == sizeof(unsigned int) )
    mpi_type_t = MPI_UNSIGNED;
  else if ( sizeof(ulong_t) == sizeof(unsigned long long) )
    mpi_type_t = MPI_UNSIGNED_LONG_LONG;
  else {
    if (comm_rank == 0) std::cout << "UNKNOWN MPI DATA TYPE." << std::endl;
    MPI_Finalize();
    return -1;
  }

  ulong_t temp_cells;
  MPI_Allreduce( &num_cells, &temp_cells, 1, mpi_type_t, MPI_SUM, MPI_COMM_WORLD );
  if ( temp_cells != global_cells ) {
    if ( comm_rank == 0) {
      std::cout << "Summed cell count should be " << global_cells
        << " got " << temp_cells << std::endl;
    }
    MPI_Finalize();
    return -1;
  }

  //============================================================================
  // Determine cell connectivity and vertex coordinates

  std::vector<real_t> coordinates(num_dims * num_verts);
  std::vector<ulong_t> cell_vertices(num_cells*num_verts_per_cell);
  std::vector<ulong_t> global_cell_id(num_cells);
  std::vector<ulong_t> global_vertex_id(num_verts);

  std::map<int,side_t> side_sets;

  //------------------------------------
  // One dimension
  if ( num_dims == 1 ) {

    auto vertex_local_id = [&](auto i) {
      return i;
    };
    auto vertex_global_id = [&](auto i) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      return istart + i;
    };

    for ( ulong_t i=0; i<local_dims[0]+1; ++i ) {
      auto id = vertex_local_id(i);
      coordinates[id              ] = local_lower[0] + i*delta[0];
      global_vertex_id[id] = vertex_global_id(i) + 1;
    }
    
    auto cell_local_id = [&](auto i) {
      return i;
    };
    auto cell_global_id = [&](auto i) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      return istart + i;
    };

    auto cell_is_left_boundary = [&](auto i)
    { return block_displ[0][block_ijk[0]] + i == 0; };
    auto cell_is_right_boundary = [&](auto i)
    { return block_displ[0][block_ijk[0]] + i == global_dims[0]-1; };

    for ( ulong_t i=0, id=0; i<local_dims[0]; ++i ) {
      auto local_id = cell_local_id(i);
      global_cell_id[local_id] = cell_global_id(i) + 1;
      cell_vertices[id++] = vertex_local_id(i  ) + 1;
      cell_vertices[id++] = vertex_local_id(i+1) + 1;
      if (cell_is_left_boundary(i)) {
        side_sets[1].name = "left";
        side_sets[1].elem_list.emplace_back(local_id+1);
        side_sets[1].side_list.emplace_back(1);
      }
      if (cell_is_right_boundary(i)) {
        side_sets[2].name = "right";
        side_sets[2].elem_list.emplace_back(local_id+1);
        side_sets[2].side_list.emplace_back(2);
      }
    }

  }
  //------------------------------------
  // Two dimensions
  else if ( num_dims == 2 ) {
    
    auto vertex_local_id = [&](auto i, auto j) {
      return i + (local_dims[0]+1)*j;
    };
    auto vertex_global_id = [&](auto i, auto j) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      auto iglobal = istart + i;
      
      auto jblk = block_ijk[1];
      auto jstart = block_displ[1][ jblk ];
      auto jglobal = jstart + j;
      return iglobal + (global_dims[0]+1)*jglobal;
    };

    for ( ulong_t j=0; j<local_dims[1]+1; ++j ) {
      for ( ulong_t i=0; i<local_dims[0]+1; ++i ) {
        auto id = vertex_local_id(i, j);
        coordinates[id              ] = local_lower[0] + i*delta[0];
        coordinates[id +   num_verts] = local_lower[1] + j*delta[1];
        global_vertex_id[id] = vertex_global_id(i,j) + 1;
      }
    }
    
    auto cell_local_id = [&](auto i, auto j) {
      return i + local_dims[0]*j;
    };
    auto cell_global_id = [&](auto i, auto j) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      auto iglobal = istart + i;
      
      auto jblk = block_ijk[1];
      auto jstart = block_displ[1][ jblk ];
      auto jglobal = jstart + j;
      return iglobal + global_dims[0]*jglobal;
    };

    auto cell_is_bottom_boundary = [&](auto, auto j)
    { return block_displ[1][block_ijk[1]] + j == 0; };
    auto cell_is_top_boundary = [&](auto, auto j)
    { return block_displ[1][block_ijk[1]] + j == global_dims[1]-1; };

    auto cell_is_left_boundary = [&](auto i, auto)
    { return block_displ[0][block_ijk[0]] + i == 0; };
    auto cell_is_right_boundary = [&](auto i, auto)
    { return block_displ[0][block_ijk[0]] + i == global_dims[0]-1; };

    for ( ulong_t j=0, id=0; j<local_dims[1]; ++j ) {
      for ( ulong_t i=0; i<local_dims[0]; ++i ) {
        auto local_id = cell_local_id(i, j);
        global_cell_id[local_id] = cell_global_id(i, j) + 1;
        cell_vertices[id++] = vertex_local_id(i  ,   j) + 1;
        cell_vertices[id++] = vertex_local_id(i+1,   j) + 1;
        cell_vertices[id++] = vertex_local_id(i+1, j+1) + 1;
        cell_vertices[id++] = vertex_local_id(i  , j+1) + 1;
        if (cell_is_bottom_boundary(i, j)) {
          side_sets[3].name = "bottom";
          side_sets[3].elem_list.emplace_back(local_id+1);
          side_sets[3].side_list.emplace_back(1);
        }
        if (cell_is_right_boundary(i, j)) {
          side_sets[2].name = "right";
          side_sets[2].elem_list.emplace_back(local_id+1);
          side_sets[2].side_list.emplace_back(2);
        }
        if (cell_is_top_boundary(i, j)) {
          side_sets[4].name = "top";
          side_sets[4].elem_list.emplace_back(local_id+1);
          side_sets[4].side_list.emplace_back(3);
        }
        if (cell_is_left_boundary(i, j)) {
          side_sets[1].name = "left";
          side_sets[1].elem_list.emplace_back(local_id+1);
          side_sets[1].side_list.emplace_back(4);
        }
      }
    }

  } // two-dimensional
  //------------------------------------
  // Three dimensions
  else if ( num_dims == 3 ) {
    
    auto vertex_local_id = [&](auto i, auto j, auto k) {
      return i + (local_dims[0]+1)*(j + (local_dims[1]+1)*k);
    };
    auto vertex_global_id = [&](auto i, auto j, auto k) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      auto iglobal = istart + i;
      
      auto jblk = block_ijk[1];
      auto jstart = block_displ[1][ jblk ];
      auto jglobal = jstart + j;
      
      auto kblk = block_ijk[2];
      auto kstart = block_displ[2][ kblk ];
      auto kglobal = kstart + k;
      return iglobal + (global_dims[0]+1)*(jglobal + (global_dims[1]+1)*kglobal);
    };
    
    auto cell_local_id = [&](auto i, auto j, auto k) {
      return i + local_dims[0]*(j + local_dims[1]*k);
    };
    auto cell_global_id = [&](auto i, auto j, auto k) {
      auto iblk = block_ijk[0];
      auto istart = block_displ[0][ iblk ];
      auto iglobal = istart + i;
      
      auto jblk = block_ijk[1];
      auto jstart = block_displ[1][ jblk ];
      auto jglobal = jstart + j;
      
      auto kblk = block_ijk[2];
      auto kstart = block_displ[2][ kblk ];
      auto kglobal = kstart + k;
      return iglobal + global_dims[0]*(jglobal + global_dims[1]*kglobal);
    };
    
    auto cell_is_left_boundary = [&](auto i, auto, auto)
    { return block_displ[0][block_ijk[0]] + i == 0; };
    auto cell_is_right_boundary = [&](auto i, auto, auto)
    { return block_displ[0][block_ijk[0]] + i == global_dims[0]-1; };
    
    auto cell_is_bottom_boundary = [&](auto, auto j, auto)
    { return block_displ[1][block_ijk[1]] + j == 0; };
    auto cell_is_top_boundary = [&](auto, auto j, auto)
    { return block_displ[1][block_ijk[1]] + j == global_dims[1]-1; };

    auto cell_is_front_boundary = [&](auto, auto, auto k)
    { return block_displ[2][block_ijk[2]] + k == 0; };
    auto cell_is_back_boundary = [&](auto, auto, auto k)
    { return block_displ[2][block_ijk[2]] + k == global_dims[2]-1; };

    for ( ulong_t k=0; k<local_dims[2]+1; ++k ) {
      for ( ulong_t j=0; j<local_dims[1]+1; ++j ) {
        for ( ulong_t i=0; i<local_dims[0]+1; ++i ) {
          auto id = vertex_local_id(i, j, k);
          coordinates[id              ] = local_lower[0] + i*delta[0];
          coordinates[id +   num_verts] = local_lower[1] + j*delta[1];
          coordinates[id + 2*num_verts] = local_lower[2] + k*delta[2];
          global_vertex_id[id] = vertex_global_id(i,j,k) + 1;
        }
      }
    }

    for ( ulong_t k=0, id=0; k<local_dims[2]; ++k ) {
      for ( ulong_t j=0; j<local_dims[1]; ++j ) {
        for ( ulong_t i=0; i<local_dims[0]; ++i ) {
          auto local_id = cell_local_id(i, j, k);
          global_cell_id[local_id] = cell_global_id(i, j, k) + 1;
          cell_vertices[id++] = vertex_local_id(i  ,   j, k) + 1;
          cell_vertices[id++] = vertex_local_id(i+1,   j, k) + 1;
          cell_vertices[id++] = vertex_local_id(i+1, j+1, k) + 1;
          cell_vertices[id++] = vertex_local_id(i  , j+1, k) + 1;

          cell_vertices[id++] = vertex_local_id(i  ,   j, k+1) + 1;
          cell_vertices[id++] = vertex_local_id(i+1,   j, k+1) + 1;
          cell_vertices[id++] = vertex_local_id(i+1, j+1, k+1) + 1;
          cell_vertices[id++] = vertex_local_id(i  , j+1, k+1) + 1;
          
          if (cell_is_bottom_boundary(i, j, k)) {
            side_sets[3].name = "bottom";
            side_sets[3].elem_list.emplace_back(local_id+1);
            side_sets[3].side_list.emplace_back(1);
          }
          if (cell_is_right_boundary(i, j, k)) {
            side_sets[2].name = "right";
            side_sets[2].elem_list.emplace_back(local_id+1);
            side_sets[2].side_list.emplace_back(2);
          }
          if (cell_is_top_boundary(i, j, k)) {
            side_sets[4].name = "top";
            side_sets[4].elem_list.emplace_back(local_id+1);
            side_sets[4].side_list.emplace_back(3);
          }
          if (cell_is_left_boundary(i, j, k)) {
            side_sets[1].name = "left";
            side_sets[1].elem_list.emplace_back(local_id+1);
            side_sets[1].side_list.emplace_back(4);
          }
          if (cell_is_front_boundary(i, j, k)) {
            side_sets[5].name = "south";
            side_sets[5].elem_list.emplace_back(local_id+1);
            side_sets[5].side_list.emplace_back(5);
          }
          if (cell_is_back_boundary(i, j, k)) {
            side_sets[6].name = "north";
            side_sets[6].elem_list.emplace_back(local_id+1);
            side_sets[6].side_list.emplace_back(6);
          }
        }
      }
    }
  
  } // three-dimensional

  // check if we need 64 bit integers
  constexpr auto max_int = std::numeric_limits<int>::max();
  auto is_int64 = 
    force_int64 ||
    (cell_vertices.size() >= max_int/2) ||
    global_verts >= max_int ||
    global_cells >= max_int;

  // size of floating point variables used in app.
  int app_word_size = sizeof(real_t);

  // size of floating point to be stored in file.
  // change to float to save space
  int exo_word_size = sizeof(real_t);

  // determine the file creation mode
  auto cmode = EX_CLOBBER;
  if (is_int64) {
    cmode |= EX_ALL_INT64_DB | EX_ALL_INT64_API;
    if (comm_rank==0) std::cout << "Using 64-bit integers." << std::endl;
  }
  else
    if (comm_rank==0) std::cout << "Using 32-bit integers." << std::endl;

  // figure out mesh file name
  std::string filename;
  if (comm_size == 1 )
    filename = name;
  else {
    auto number = comm_size;
    unsigned int num_digits = 1;
    while (number) {
      number /= 10;
      num_digits++;
    }
    std::stringstream ss;
    ss
      << name
      << "."
      << std::setfill('0') << std::setw(num_digits) << comm_size
      << "."
      << std::setfill('0') << std::setw(num_digits) << comm_rank;
    filename = ss.str();
  }

  // create file
  auto exo_id =
      ex_create(filename.c_str(), cmode, &app_word_size, &exo_word_size);
  if (exo_id < 0) {
    if ( comm_rank == 0 ) {
      std::cerr << "Problem writing exodus file, ex_create() returned " << exo_id << std::endl;
    }
    MPI_Finalize();
    return exo_id;
  }
  else {
    if (comm_rank==0) std::cout << "Opened file for writing: " << name << std::endl;
  }

  
  ex_init_params exopar;
  std::strcpy(exopar.title, "Exodus II output from exogen.");
  exopar.num_dim = num_dims;
  exopar.num_nodes = num_verts;
  exopar.num_edge = 0;
  exopar.num_edge_blk = 0;
  exopar.num_face = 0;
  exopar.num_face_blk = 0;
  exopar.num_elem = num_cells;
  exopar.num_elem_blk = 1;
  exopar.num_node_sets = 0;
  exopar.num_edge_sets = 0;
  exopar.num_face_sets = 0;
  exopar.num_side_sets = side_sets.size();
  exopar.num_elem_sets = 0;
  exopar.num_node_maps = 0;
  exopar.num_edge_maps = 0;
  exopar.num_face_maps = 0;
  exopar.num_elem_maps = 0;

  auto status = ex_put_init_ext(exo_id, &exopar);
  if (status) {
    if (comm_rank==0)
      std::cerr << "Problem putting exodus file parameters, ex_put_init_ext() returned "
          << status << std::endl;
    MPI_Finalize();
    return status;
  }

  // exodus is kind enough to fetch the data in the real type we ask for
  status = ex_put_coord(
      exo_id, coordinates.data(), coordinates.data() + num_verts,
      coordinates.data() + 2 * num_verts);

  if (status) {
    if (comm_rank==0)
      std::cerr
          << "Problem putting vertex coordinates to exodus file, "
          << " ex_put_coord() returned " << status << std::endl;
    MPI_Finalize();
    return status;
  }
  else {
    if (comm_rank==0)
      std::cout << "Wrote " << global_verts << " vertex coordinates" << std::endl;
  }

  const char *coord_names[3] = {"x", "y", "z"};
  status = ex_put_coord_names (exo_id, const_cast<char**>(coord_names));
  if (status) {
    if (comm_rank==0)
      std::cerr
          << "Problem putting coordinate names to exodus file, "
          << " ex_put_coord_names() returned " << status << std::endl;
    MPI_Finalize();
    return status;
  }

  ex_entity_id elem_blk_id = 1;
  if (is_int64) 
    status = write_cells<long long>(exo_id, elem_blk_id, cell_vertices, num_verts_per_cell);
  else
    status = write_cells<int>(exo_id, elem_blk_id, cell_vertices, num_verts_per_cell);
  if (status) {
    if (comm_rank==0)
      std::cerr
          << "Problem putting elements to exodus file, "
          << " ex_put_block() returned " << status << std::endl;
    MPI_Finalize();
    return status;
  }
  else {
    if (comm_rank==0) std::cout << "Wrote " << global_cells << " cells" << std::endl;
  }

  if (is_int64) 
    status = write_mapping<long long>(exo_id, global_cell_id, global_vertex_id);
  else
    status = write_mapping<int>(exo_id, global_cell_id, global_vertex_id);
  if (status) {
    if (comm_rank==0)
      std::cerr
          << "Problem putting mapping to exodus file, "
          << " ex_put_num_map() returned " << status << std::endl;
    MPI_Finalize();
    return status;
  }
  else {
    if (comm_rank==0) std::cout << "Wrote mapping." << std::endl;
  }
 
  std::vector<char *> ss_names;
  ss_names.reserve(side_sets.size());

  for (const auto & ss : side_sets) {
    ss_names.emplace_back( const_cast<char*>(ss.second.name.c_str()) );
    const auto & elem_list = ss.second.elem_list;
    const auto & side_list = ss.second.side_list;
    auto sid = ss.first;
    if (is_int64)
      status = write_side_set<long long> (exo_id, sid, elem_list, side_list);
    else
      status = write_side_set<int> (exo_id, sid, elem_list, side_list);
    if (status) {
      if (comm_rank==0)
        std::cerr
            << "Problem putting side sets to exodus file, "
            << " ex_put_side_set() returned " << status << std::endl;
      MPI_Finalize();
      return status;
    }
    sid++;
  }

  status = ex_put_names(exo_id, EX_SIDE_SET, ss_names.data());
  if (status) {
    if (comm_rank==0)
      std::cerr
        << "Problem putting side set names, ex_put_names() returned " << exo_id
        << std::endl;
    MPI_Finalize();
    return status;
  }

    
  status = ex_close(exo_id);
  if (status) {
    if (comm_rank==0)
      std::cerr
        << "Problem closing exodus file, ex_close() returned " << exo_id
        << std::endl;
    MPI_Finalize();
    return status;
  }
  else {
    if (comm_rank==0) std::cout << "Closed exodus file!" << std::endl;
  }

  return MPI_Finalize();

}
