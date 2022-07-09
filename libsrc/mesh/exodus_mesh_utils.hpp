#ifndef EXODUS_MESH_UTILS_HPP
#define EXODUS_MESH_UTILS_HPP

#include "comm/comm_utils.hpp"
#include "comm/mpi_comm.hpp"
#include "io/exo_reader.hpp"
#include "io/side_set.hpp"
#include "math/subdivide.hpp"
#include "utils/crs.hpp"

#define CHUNK_SIZE 256

namespace prl {

//==============================================================================
// Lambda function to process a block
//==============================================================================
template< typename T, typename U, typename V, typename W>
static void filter_block(
    const T & counts_in,
    const U & indices_in,
    size_t min,
    size_t max,
    size_t & counter,
    V & offsets_out,
    W & indices_out)
{
  auto num_offsets = counts_in.size();
  auto num_indices = indices_in.size();

  // if this is outside the range, just increment counter and return
  if ( counter + num_offsets <= min || counter > max )
  {
    counter += num_offsets;
    return;
  }

  // storage for element verts
  offsets_out.reserve( offsets_out.size() + num_offsets );
  indices_out.reserve( indices_out.size() + num_indices );

  if (offsets_out.empty()) offsets_out.emplace_back(0);
  
  // create cells in mesh
  size_t base = 0;
  for (size_t e = 0; e < num_offsets; ++e, ++counter) {
    // get the number of nodes
    auto cnt = counts_in[e];
    // the global id is within the range
    if ( counter >= min && counter <= max )
    {
      // copy local vertices into vector ( exodus uses 1 indexed arrays )
      for (int v = 0; v < cnt; v++)
        indices_out.emplace_back(indices_in[base + v] - 1);
      // add the row
      offsets_out.emplace_back( offsets_out.back() + cnt );
    }
    // base offset into elt_conn
    base += cnt;
  }
}

//============================================================================
// filter out sides
//============================================================================
template< typename T >
void filter_sides(
    size_t ss_id,
    const std::vector<T> & side_set_elem_list,
    const std::vector<T> & side_set_side_list,
    size_t cell_min,
    size_t cell_max,
    std::vector<int_t> & side_id,
    std::vector<int_t> & side_cell,
    std::vector<int_t> & side_cell_face )
{

  auto num_side_in_set = side_set_elem_list.size();
  
  // filter sides for elementes i own
  std::vector<size_t> vs;
  for ( size_t j=0; j<num_side_in_set; ++j ) {
    size_t global_id = side_set_elem_list[j] - 1;
    if (global_id >= cell_min && global_id <= cell_max )
    {
      // add side info
      side_id.emplace_back( ss_id-1 );
      side_cell.emplace_back( global_id );
      side_cell_face.emplace_back( side_set_side_list[j] - 1 );
    }
  }

}
  
  
//============================================================================
// Filter out sides
//============================================================================
template< typename T >
static void filter_sides(
    size_t ss_id,
    const std::vector<T> & side_set_elem_list,
    const std::vector<T> & side_set_side_list,
    size_t cell_min,
    size_t cell_max,
    const std::map<long_t, int_t> & cell_global2local,
    std::vector<int_t> & side_id,
    std::vector<int_t> & side_cell,
    std::vector<int_t> & side_cell_face )
{

  auto num_side_in_set = side_set_elem_list.size();

  // filter sides for elementes i own
  for ( size_t j=0; j<num_side_in_set; ++j ) {
    size_t elem_gid = side_set_elem_list[j] - 1;
    if (elem_gid >= cell_min && elem_gid <= cell_max )
    {
      auto elem_lid = cell_global2local.at(elem_gid);
      // add side info
      side_id.emplace_back( ss_id-1 );
      side_cell.emplace_back( elem_lid );
      side_cell_face.emplace_back( side_set_side_list[j]-1 );
    }
  }

}
//============================================================================
//! \brief read cell ids
//============================================================================
static int read_exo_cell_ids(
    exo_reader_t & exo,
    bool serial,
    bool int64,
    size_t cell_min,
    size_t num_cells,
    std::vector<long_t> & cell_local2global,
    std::map<long_t, int_t> & cell_global2local)
{
  
  int ret = 0;

  // for serial files, make up the indices
  if ( serial ) {
    cell_local2global.clear();
    cell_local2global.resize(num_cells);
    for ( size_t i=0; i<num_cells; ++i ) {
      cell_local2global[i] = cell_min + i;
    }
  }
  // read them for parallel files
  else if (int64) {
    ret = exo.read_element_map<long long>(num_cells, cell_local2global);
  }
  else {
    ret = exo.read_element_map<int>(num_cells, cell_local2global);
  }

  if (ret) return ret;
  
  // invert the global id to local id map
  for ( size_t i=0; i<num_cells; ++i ) {
    auto global_id = cell_local2global[i];
    cell_global2local.emplace( std::make_pair( global_id, i ) );
  }

  return 0;

}

//============================================================================
//! \brief read a block
//============================================================================
template <typename U>
int read_exo_block(
    exo_reader_t & exo,
    ex_entity_id blk_id,
    size_t cell_min,
    size_t cell_max,
    const std::vector<size_t> & cell_partitioning,
    mpi_comm_t & comm,
    bool do_read,
    bool do_send,
    size_t & cell_counter,
    crs_t<long_t, int_t> & cell2vertices,
    shape_t & block_type)
{

  auto comm_rank = comm.rank(); 
  auto comm_size = comm.size();

  // get the info about this block
  using stats_t = exo_block_stats_t<U>;
  stats_t block_stats;
  if  (do_read) {
    auto stat = exo.read_block_stats(
      blk_id,
      EX_ELEM_BLOCK,
      block_stats);
    if (stat) {
      std::cout << "Couldn't read block stats." << std::endl;
      return stat;
    }
  }
  if (do_send) comm.broadcast(block_stats, 0);
  
  if (!block_stats.num_elem_this_blk) {
    block_type = shape_t::empty;
    return 0;
  }

  
  //----------------------------------
  // Serial version
  if (do_read && !do_send) {
  
    std::vector<int> counts;
    std::vector<U> indices;
  
    bool finished;
    auto res = exo.read_block<U>(
        blk_id,
        EX_ELEM_BLOCK,
        block_stats,
        counts,
        indices,
        block_type,
        finished);
    if (res) {
      std::cout << "Couldn't read block." << std::endl;
      return res;
    }
    filter_block(
        counts,
        indices,
        cell_min,
        cell_max,
        cell_counter,
        cell2vertices.offsets,
        cell2vertices.indices );
    return 0;
  }

  //----------------------------------
  // Parallel version
  else {
  
    block_type = shape_t::empty;

    struct send_data_t {
      shape_t block_type;
      size_t num_counts;
      size_t num_indices;
    };
    
    std::vector<int> counts;
    std::vector<U> indices;
    send_data_t data;
      
    mpi_request_t requests;

    double elapsed_time = 0;
    size_t elapsed_counter = 0;

    size_t first_read_cell = cell_counter;
    size_t last_read_cell = cell_counter + block_stats.num_elem_this_blk;

    std::vector<size_t> distribution(comm_size+1);
    distribution[0] = 0;

    for (int r=0; r<comm_size; ++r) {
      auto begin = std::max(first_read_cell, cell_partitioning[r]);
      auto end = std::min(last_read_cell, cell_partitioning[r+1]);
      size_t sz = end>begin ? end-begin : 0;
      distribution[r+1] = distribution[r] + sz;
    }

    for (int r=1; r<comm_size; ++r) {
      size_t start = distribution[r];
      auto rank_cells = distribution[r+1] - distribution[r];
    
      while (rank_cells > 0) {

        auto chunk = std::min<size_t>(CHUNK_SIZE, rank_cells);

        //----------------------------
        // Root
        if (comm_rank == 0) {
          
          bool finished;
          auto status = exo.read_block<U>(
              blk_id,
              EX_ELEM_BLOCK,
              block_stats,
              counts,
              indices,
              block_type,
              finished,
              start,
              chunk);
          if (status) return status;

          data = {block_type, counts.size(), indices.size()};
          
          auto req1 = comm.send(data, r, 0);
          auto req2 = comm.send(counts, r, 1);
          auto req3 = comm.send(indices, r, 2);
          requests.transfer(req1);
          requests.transfer(req2);
          requests.transfer(req3);
              
          auto startt = MPI_Wtime();
          requests.wait();
          elapsed_time += MPI_Wtime() - startt;
          elapsed_counter++;
        }
        //----------------------------
        // Target
        else if (comm_rank == r)  {
          auto counts_size = counts.size();
          auto indices_size = indices.size();

          auto req1 = comm.receive(data, 0, 0);
          req1.wait();

          block_type = data.block_type;
          counts.resize(counts_size + data.num_counts);
          indices.resize(indices_size + data.num_indices);
          
          auto counts_ref = make_array_ref(
              counts.data() + counts_size,
              data.num_counts);
          auto req2 = comm.receive(counts_ref, 0, 1);

          auto indices_ref = make_array_ref(
              indices.data() + indices_size,
              data.num_indices);
          auto req3 = comm.receive(indices_ref, 0, 2);
            
          requests.transfer(req2);
          requests.transfer(req3);
          requests.wait();
    
        }
        // End  root / target
        //----------------------------

        start += chunk;
        rank_cells -= chunk;
      } // while
      
      // without this, some ranks seem to get too far ahead and timeout
      comm.barrier();    

    } // ranks
      
    if (comm_rank == 0 && elapsed_counter) {
      std::cout << "Block<" << blk_id << ">: Spent a total/average of " << elapsed_time
        << " / " << elapsed_time / elapsed_counter << " secs waiting." << std::endl;
    }
      
    // root needs to read its data now
    if (comm_rank ==0 ) {
      size_t start = distribution[0];
      auto rank_cells = distribution[1] - distribution[0];
      indices.clear();
      counts.clear();
      if (rank_cells) {
        bool finished;
        auto res = exo.read_block<U>(
            blk_id,
            EX_ELEM_BLOCK,
            block_stats,
            counts,
            indices,
            block_type,
            finished,
            start,
            rank_cells);
        if (res) return res;
      }
    }

    // now filter cells
    auto rank_size = distribution[comm_rank+1] - distribution[comm_rank];
    if (rank_size) {
      cell_counter += distribution[comm_rank];
      filter_block(
          counts,
          indices,
          cell_min,
          cell_max,
          cell_counter,
          cell2vertices.offsets,
          cell2vertices.indices );
    }

    // make sure everyone has the last cell counter
    cell_counter = last_read_cell;
    
    return 0;
  }
}

//============================================================================
//! \brief read the coordinates of the mesh from a file.
//! \param [in] exo_id  The exodus file id.
//! \return the vertex coordinates
//============================================================================
static int read_exo_point_coords_serial(
    exo_reader_t & exo,
    const ex_init_params & exo_params,
    mpi_comm_t & comm,
    bool do_send,
    const crs_t<long_t, int_t> & cell2vertices,
    std::map<long_t, int_t> & vertex_global2local,
    std::vector<long_t> & vertex_local2global,
    std::vector<real_t> & vertices)
{
  
  auto comm_rank = comm.rank(); 
  auto comm_size = comm.size();

  // create a local numbering of the vertices
  size_t local_vertices{0};
  vertex_global2local.clear();
 
  size_t num_cells = cell2vertices.size();
  for ( size_t i=0; i<num_cells; ++i ) {
    for ( auto j=cell2vertices.offsets[i]; j<cell2vertices.offsets[i+1]; ++j ) {
      auto global_id = cell2vertices.indices[j];
      auto res = vertex_global2local.emplace(
          std::make_pair(global_id, local_vertices) );
      if ( res.second ) local_vertices++;
    }
  }
  
  size_t num_dims = exo_params.num_dim;
  size_t num_nodes = exo_params.num_nodes;

  vertices.clear();
  vertices.resize(local_vertices*num_dims);

  //--------------------------------
  // Read vertices in serial
  if (!do_send) {
    std::vector<real_t> coordinates;
    auto ret = exo.read_point_coords(num_dims, num_nodes, coordinates);
    if (ret) return ret;

    for ( const auto & id_pair : vertex_global2local ) {
      auto global_id = id_pair.first;
      auto local_id = id_pair.second;
      for ( unsigned d=0; d<num_dims; ++d ) {
        vertices[ local_id*num_dims + d ] = 
          coordinates[d * local_vertices + global_id];
      }
    }

    // invert the global id to local id map
    vertex_local2global.clear();
    vertex_local2global.resize(local_vertices);
      
    for ( const auto & global_local : vertex_global2local )
      vertex_local2global[global_local.second] = global_local.first;

  }
  //--------------------------------
  // Read vertices in chunks
  else {
    
    std::vector<size_t> distribution;
    subdivide(num_nodes, comm_size, distribution);

    std::vector<real_t> coordinates, vert_buf;
  
    double elapsed_time = 0;
    size_t elapsed_counter = 0;

    auto comm_verts = distribution[comm_rank+1] - distribution[comm_rank];
    if (comm_rank != 0) coordinates.resize(comm_verts*num_dims);

    for (int r=1; r<comm_size; ++r) {
      
      size_t start = distribution[r];
      auto rank_verts = distribution[r+1] - distribution[r];

      while (rank_verts > 0) {

        auto chunk = std::min<size_t>(CHUNK_SIZE, rank_verts);
        
        //----------------------------
        // Root
        if (comm_rank == 0) {

          auto ret = exo.read_point_coords(
              num_dims, start, start+chunk, coordinates);
          if (ret) return ret;
      
          auto req = comm.send(coordinates, r, 0);
          
          auto startt = MPI_Wtime();
          req.wait();
          elapsed_time += MPI_Wtime() - startt;
          elapsed_counter++;
        }
        
        //----------------------------
        // Target
        else if (comm_rank == r)  {
          
          vert_buf.resize(chunk*num_dims);
          comm.receive(vert_buf, 0, 0);

          auto offset = start - distribution[r];
          for (size_t i=0; i<chunk; ++i) {
            auto pos = offset + i;
            for (unsigned d=0; d<num_dims; ++d) {
              coordinates[d*comm_verts + pos] = vert_buf[d*chunk + i];
            }
          }
        
        }
        // End  root / target
        //----------------------------

        start += chunk;
        rank_verts -= chunk;
      } // while
      
      // without this, some ranks seem to get too far ahead and timeout
      comm.barrier();
        
    } // ranks
    
    if (comm_rank == 0 && elapsed_counter) {
      std::cout << "Vertices: Spent a total/average of " << elapsed_time 
        << " / " << elapsed_time / elapsed_counter << " secs waiting." << std::endl;
    }

    // root needs to read its data
    if (comm_rank == 0) {
      size_t start = distribution[0];
      auto rank_verts = distribution[1] - distribution[0];
      if (rank_verts) {
        auto ret = exo.read_point_coords(
            num_dims,
            start,
            start+rank_verts,
            coordinates);
        if (ret) return ret;
      }
    }

    // invert the global id to local id map
    vertex_local2global.clear();
    vertex_local2global.resize(local_vertices);
    for ( const auto & global_local : vertex_global2local )
      vertex_local2global[global_local.second] = global_local.first;
    
    // figure out who owns the vertices i need
    std::vector<int> rank_owner(local_vertices);
    
    // determine the rank owners and send counts
    std::vector<int> sendcounts(comm_size, 0);

    int r = 0;
    for ( const auto & global_local : vertex_global2local ) {
      size_t global_id = global_local.first;
      auto local_id = global_local.second;
      while (global_id >= distribution[r+1]) { ++r; }
      rank_owner[local_id] = r;
      sendcounts[r]++;
    }

    std::vector<int> recvcounts(comm_size);
    comm.all_to_all(sendcounts, recvcounts);
  
    // finish displacements
    std::vector<int> senddispls(comm_size+1);
    std::vector<int> recvdispls(comm_size+1);
    senddispls[0] = 0;
    recvdispls[0] = 0;
    for(int r = 0; r < comm_size; ++r)  {
      senddispls[r + 1] = senddispls[r] + sendcounts[r];
      recvdispls[r + 1] = recvdispls[r] + recvcounts[r];
    }

    std::vector<size_t> index_buf(senddispls[comm_size]);
    std::fill(sendcounts.begin(), sendcounts.end(), 0);

    for ( const auto & global_local : vertex_global2local ) {
      auto global_id = global_local.first;
      auto local_id = global_local.second;
      auto r = rank_owner[local_id];
      auto pos = senddispls[r] + sendcounts[r];
      index_buf[pos] = global_id;
      sendcounts[r]++;
    }

    // now send the vertices i need to the owners
    std::vector<size_t> send_indices(recvdispls[comm_size]);
    comm.all_to_allv(
        index_buf,
        sendcounts,
        senddispls,
        send_indices,
        recvcounts,
        recvdispls);
  
    // copy ther vertices
    auto start = distribution[comm_rank];
    vert_buf.resize(send_indices.size()*num_dims);

    for (size_t i=0; i<send_indices.size(); ++i) {
      auto pos = send_indices[i] - start;
      for (unsigned d=0; d<num_dims; ++d)
        vert_buf[i*num_dims + d] = coordinates[d*comm_verts + pos];
    }

    // swap send and recv counts
    std::swap(sendcounts, recvcounts);
    std::swap(senddispls, recvdispls);

    for (int i=0; i<comm_size; ++i) {
      sendcounts[i] *= num_dims;
      senddispls[i+1] *= num_dims;
      recvcounts[i] *= num_dims;
      recvdispls[i+1] *= num_dims;
    }
    
    coordinates.resize(local_vertices*num_dims);
    comm.all_to_allv(
        vert_buf,
        sendcounts,
        senddispls,
        coordinates,
        recvcounts,
        recvdispls);
    
    // unpack
    for (int i=0; i<comm_size; ++i) recvdispls[i+1] /= num_dims;
    
    std::fill(recvcounts.begin(), recvcounts.end(), 0);

    for ( const auto & global_local : vertex_global2local ) {
      //auto global_id = global_local.first;
      auto local_id = global_local.second;
      auto r = rank_owner[local_id];
      auto pos = recvdispls[r] + recvcounts[r];
      // here
      for ( unsigned d=0; d<num_dims; ++d )
        vertices[ local_id*num_dims + d ] = coordinates[ pos*num_dims + d ];
      recvcounts[r]++;
    }

  }

  return 0;
    
}
  
//============================================================================
//! \brief read the coordinates of the mesh from a file.
//! \param [in] exo_id  The exodus file id.
//! \return the vertex coordinates
//============================================================================
static int read_exo_point_coords_parallel(
    exo_reader_t & exo,
    const ex_init_params & exo_params,
    bool int64,
    std::map<long_t, int_t> & vertex_global2local,
    std::vector<long_t> & vertex_local2global,
    std::vector<real_t> & vertices)
{

  std::size_t num_dims = exo_params.num_dim;
  std::size_t num_nodes = exo_params.num_nodes;
  
  std::vector<real_t> coords;
  auto ret = exo.read_point_coords(
      num_dims,
      num_nodes,
      coords);
  if (ret) return ret;
  if (coords.size() != num_dims * num_nodes)
    return -1;

  auto num_vertices = coords.size() / num_dims;
 
  if (int64) {
    ret = exo.read_node_map<long long>(num_vertices, vertex_local2global);
  }
  else {
    ret = exo.read_node_map<int>(num_vertices, vertex_local2global);
  }
  if (ret) return ret;
  
  vertices.clear();
  vertices.resize(num_vertices*num_dims);

  for ( size_t i=0; i<num_vertices; ++i ) {
    auto global_id = vertex_local2global[i];
    vertex_global2local.emplace( std::make_pair( global_id, i ) );
    for (unsigned d=0; d<num_dims; ++d)
      vertices[i*num_dims + d] = coords[d*num_vertices + i];
  }

  return 0;

} 

//============================================================================
//! \brief read the faces
//============================================================================
template<typename U>
int read_exo_face_block(
    exo_reader_t & exo,
    ex_entity_id blk_id,
    mpi_comm_t & comm,
    bool do_read,
    bool do_send,
    const crs_t<long_t, int_t> & cell2entities,
    int_t num_regions,
    const std::vector<int_t> & region_cell_offsets,
    const std::vector<shape_t> & cell_type,
    size_t face_start,
    size_t & num_faces,
    std::map<long_t, int_t> & face_global2local,
    crs_t<long_t, int_t> & face2vertices,
    std::map<std::vector<long_t>, long_t> sorted_face2vertices)
{

  int ret;
  auto comm_rank = comm.rank(); 
  auto comm_size = comm.size();

  //----------------------------------------------------------------------------
  // Read face block stats

  using stats_t = exo_block_stats_t<U>;
  stats_t block_stats;
    
  // get the block statistics
  if (do_read) {
    ret = exo.read_block_stats(
        blk_id,
        EX_FACE_BLOCK,
        block_stats);
    if (ret) {
      std::cout << "Couldn't read face block stats." << std::endl;
      return ret;
    }
  }
  if (do_send) comm.broadcast(block_stats, 0);
    
  num_faces = block_stats.num_elem_this_blk;
  auto face_end = face_start + num_faces;
    
  // read each block (add all faces right now.  we will
  // filter out the unused ones later
  bool finished;
  std::vector<int> counts;
  std::vector<U> indices;
  shape_t block_type;
 

  //----------------------------------------------------------------------------
  // Read vertices in serial
  if (do_read && !do_send) {
  
    // first get the vertices
    ret = exo.read_block<U>(
        blk_id,
        EX_FACE_BLOCK,
        block_stats,
        counts,
        indices,
        block_type,
        finished);
    if (ret) {
      std::cout << "Couldn't read face block." << std::endl;
      return ret;
    }
    
    // Subtract one fromn vertices
    for (auto & i : indices) i--;
  
    // insert face vertices
    std::vector<long_t> sorted_vs;
    
    for (size_t i=0, j=0; i<counts.size(); ++i) {
      auto global_id = i + face_start;
      auto n = counts[i];
      // sort vertices
      sorted_vs.assign(&indices[j], &indices[j+n]);
      std::sort(sorted_vs.begin(), sorted_vs.end());
      // try adding face
      auto local_id = sorted_face2vertices.size();
      auto res = sorted_face2vertices.emplace( sorted_vs, local_id );
      // now add mapping if new
      if (res.second) face2vertices.push_back(&indices[j], &indices[j+n]);
      face_global2local.emplace(global_id, res.first->second);
      j += n;
    }

  }
  //----------------------------------------------------------------------------
  // Read vertices in chunks
  else {
    
    std::vector<size_t> distribution;
    subdivide(num_faces, comm_size, distribution);
    auto rank_faces = distribution[comm_rank+1] - distribution[comm_rank];

    counts.reserve(rank_faces+1);

    struct send_data_t {
      shape_t block_type;
      size_t num_counts;
      size_t num_indices;
    };
    
    send_data_t data;
    
    mpi_request_t requests;
  
    double elapsed_time = 0;
    size_t elapsed_counter = 0;

    for (int r=1; r<comm_size; ++r) {
      
      size_t start = distribution[r];
      auto faces_left = distribution[r+1] - distribution[r];

      while (faces_left > 0) {

        auto chunk = std::min<size_t>(CHUNK_SIZE, faces_left);
        
        //----------------------------
        // Root
        if (comm_rank == 0) {

          auto ret = exo.read_block<U>(
              blk_id,
              EX_FACE_BLOCK,
              block_stats,
              counts,
              indices,
              block_type,
              finished,
              start,
              chunk);
          if (ret) return ret;
      
          auto req1 = comm.send(data, r, 0);
          auto req2 = comm.send(counts, r, 1);
          auto req3 = comm.send(indices, r, 2);
          requests.transfer(req1);
          requests.transfer(req2);
          requests.transfer(req3);
              
          auto startt = MPI_Wtime();
          requests.wait();
          elapsed_time += MPI_Wtime() - startt;
          elapsed_counter++;
        }
        
        //----------------------------
        // Target
        else if (comm_rank == r)  {
          
          auto counts_size = counts.size();
          auto indices_size = indices.size();

          auto req1 = comm.receive(data, 0, 0);
          req1.wait();

          block_type = data.block_type;
          counts.resize(counts_size + data.num_counts);
          indices.resize(indices_size + data.num_indices);
          
          auto counts_ref = make_array_ref(
              counts.data() + counts_size,
              data.num_counts);
          auto req2 = comm.receive(counts_ref, 0, 1);

          auto indices_ref = make_array_ref(
              indices.data() + indices_size,
              data.num_indices);
          auto req3 = comm.receive(indices_ref, 0, 2);
            
          requests.transfer(req2);
          requests.transfer(req3);
          requests.wait();
        
        }
        // End  root / target
        //----------------------------

        start += chunk;
        faces_left -= chunk;
      } // while
      
      // without this, some ranks seem to get too far ahead and timeout
      comm.barrier();
        
    } // ranks
    
    if (comm_rank == 0 && elapsed_counter) {
      std::cout << "Face block<" << blk_id << ">: Spent a total/average of ";
      std::cout << elapsed_time << " / " << elapsed_time / elapsed_counter;
      std::cout << " secs waiting." << std::endl;
    }
      
    // root needs to read its data
    if (comm_rank == 0 && rank_faces) {
      size_t start = distribution[0];
      auto res = exo.read_block<U>(
          blk_id,
          EX_FACE_BLOCK,
          block_stats,
          counts,
          indices,
          block_type,
          finished,
          start,
          rank_faces);
      if (res) return res;
    }
    
    // subtract one from indices
    for (auto & i : indices) i--;

    // turn counts into offsets
    counts.insert(counts.begin(), 0);
    std::partial_sum(counts.begin(), counts.end(), counts.begin());
    
    //--------------------------------------------------------------------------
    // figure out who owns the faces i need
    
    // determine the rank owners and send counts
    std::vector<int> rank_owner;
    std::vector<int> sendcounts(comm_size, 0);

    for (int_t r=0; r<num_regions; ++r) {
      auto region_start = region_cell_offsets[r];
      auto region_end = region_cell_offsets[r+1];
      auto elem_type = cell_type[region_start];
      if (elem_type != shape_t::polyhedron) continue;
      for (auto c=region_start; c<region_end; ++c) {
        auto fs = cell2entities.at(c);
        for (size_t f : fs) {
          if (f >= face_start && f < face_end) {
            size_t fid = f - face_start;
            auto rank = owner(distribution, fid); 
            rank_owner.emplace_back(rank);
            sendcounts[rank]++;
          }
        } // fs
      } // cells
    } // regions
    
    std::vector<int> recvcounts(comm_size);
    comm.all_to_all(sendcounts, recvcounts);
  
    // finish displacements
    std::vector<int> senddispls(comm_size+1);
    std::vector<int> recvdispls(comm_size+1);
    senddispls[0] = 0;
    recvdispls[0] = 0;
    std::partial_sum(sendcounts.begin(), sendcounts.end(), &senddispls[1]);
    std::partial_sum(recvcounts.begin(), recvcounts.end(), &recvdispls[1]);

    std::vector<size_t> sendbuf(senddispls[comm_size]);
    std::fill(sendcounts.begin(), sendcounts.end(), 0);

    for (int_t r=0, face_counter=0; r<num_regions; ++r) {
      auto region_start = region_cell_offsets[r];
      auto region_end = region_cell_offsets[r+1];
      auto elem_type = cell_type[region_start];
      if (elem_type != shape_t::polyhedron) continue;
      for (auto c=region_start; c<region_end; ++c) {
        auto fs = cell2entities.at(c);
        for (size_t f : fs) {
          if (f >= face_start && f < face_end) {
            size_t fid = f - face_start;
            auto rank = rank_owner[face_counter];
            auto pos = senddispls[rank] + sendcounts[rank];
            sendbuf[pos] = fid;
            sendcounts[rank]++;
            face_counter++;
          }
        } // fs
      } // cells
    } // regions

    // now send the vertices i need to the owners
    std::vector<size_t> send_indices(recvdispls[comm_size]);
    comm.all_to_allv(
        sendbuf,
        sendcounts,
        senddispls,
        send_indices,
        recvcounts,
        recvdispls);
    
    //--------------------------------------------------------------------------
    // Send back the results
    
    // send the counts
    std::fill(sendcounts.begin(), sendcounts.end(), 0);
    
    for (int_t r=0; r<comm_size; ++r) {
      for (auto i=recvdispls[r]; i<recvdispls[r+1]; ++i) {
        // send id, count and indices
        auto fid = send_indices[i];
        auto n = counts[fid+1] - counts[fid];
        sendcounts[r] += 2 + n;
      }
    }
    
    comm.all_to_all(sendcounts, recvcounts);
    
    // finish dispacements
    std::partial_sum(sendcounts.begin(), sendcounts.end(), &senddispls[1]);
    std::partial_sum(recvcounts.begin(), recvcounts.end(), &recvdispls[1]);

    // now fill the buffer
    sendbuf.clear();
    sendbuf.resize( senddispls.back() );

    for (size_t i=0, cnt=0; i<send_indices.size(); ++i) {
      // send count and indices
      auto fid = send_indices[i];
      auto start = counts[fid];
      auto end = counts[fid+1];
      sendbuf[cnt] = fid;
      cnt++;
      sendbuf[cnt] = end - start;
      cnt++;
      for (int j=start; j<end; ++j) {
        sendbuf[cnt] = indices[j]; 
        cnt++;
      } // faces
    } // indices
    

    std::vector<size_t> recvbuf(recvdispls.back());
    comm.all_to_allv(
        sendbuf,
        sendcounts,
        senddispls,
        recvbuf,
        recvcounts,
        recvdispls);

    //--------------------------------------------------------------------------
    // Unpack
    
    std::vector<long_t> vs;
    for (size_t cnt=0; cnt<recvbuf.size();) {
      auto global_id = recvbuf[cnt] + face_start;
      cnt++;
      auto n = recvbuf[cnt];
      cnt++;
      vs.assign(&recvbuf[cnt], &recvbuf[cnt+n]);
      std::sort(vs.begin(), vs.end());
      auto local_id = sorted_face2vertices.size();
      auto res = sorted_face2vertices.emplace(vs, local_id);
      if (res.second) face2vertices.push_back(&recvbuf[cnt], &recvbuf[cnt+n]);
      face_global2local.emplace(global_id, res.first->second);
      cnt += n;
    } // cnt
  
  }
  

  return 0;
    
}
  
//============================================================================
//! \brief read a block
//============================================================================
template <typename U>
int read_exo_side_set(
    exo_reader_t & exo,
    mpi_comm_t & comm,
    ex_entity_id ss_id,
    size_t cell_min,
    size_t cell_max,
    bool do_read,
    bool do_send,
    const std::map<long_t, int_t> & element_global2local,
    std::vector<int_t> & side_id,
    std::vector<int_t> & side_cell,
    std::vector<int_t> & side_cell_face)
{
  auto comm_rank = comm.rank();

  // get the info about this block
  using stats_t = exo_side_set_stats_t<U>;
  stats_t side_set_stats;
  if  (do_read) {
    auto ret = exo.read_side_set_stats(ss_id, side_set_stats);
    if (ret) return ret;
  }
  if (do_send) comm.broadcast(side_set_stats, 0);
    
  if (!side_set_stats.num_side_in_set) return 0;
  
  bool finished;
  std::vector<U> side_set_elem_list;
  std::vector<U> side_set_side_list; 
    
  //----------------------------------
  // Serial version
  if (do_read && !do_send) {
  
    auto res = exo.read_side_set<U>(
        ss_id,
        side_set_stats,
        side_set_elem_list,
        side_set_side_list,
        finished);
    if (res) return res;

    filter_sides( 
        ss_id,
        side_set_elem_list,
        side_set_side_list,
        cell_min,
        cell_max,
        side_id,
        side_cell,
        side_cell_face );

  }
  //----------------------------------
  // Parallel version
  else {

    // now read side sets
    bool more_work = true;
    size_t start = 0;
    constexpr auto chunk = CHUNK_SIZE;
    
    struct send_data_t {
      bool more_work;
      size_t size;
    };
    
    send_data_t data;
    
    mpi_request_t requests;
    
    double elapsed_time = 0;
    size_t elapsed_counter = 0;
    
    while (more_work) {

      if (comm_rank == 0) {
        if (start) {
          auto startt = MPI_Wtime();
          requests.wait();
          elapsed_time += MPI_Wtime() - startt;
          elapsed_counter++;
        }
        bool finished;
        auto res = exo.read_side_set<U>(
            ss_id,
            side_set_stats,
            side_set_elem_list,
            side_set_side_list,
            finished,
            start,
            chunk);
        if (res) return res;
        more_work = !finished;
        data = {more_work, side_set_elem_list.size()};
      }
      
      
      auto req1 = comm.broadcast(data, 0);
      if (comm_rank != 0) req1.wait();
      else requests.transfer(req1);

      more_work = data.more_work;

      if (comm_rank != 0) {
        side_set_elem_list.clear();
        side_set_side_list.clear();
        side_set_elem_list.resize(data.size);
        side_set_side_list.resize(data.size);
      }

      auto req2 = comm.broadcast(side_set_elem_list, 0);
      auto req3 = comm.broadcast(side_set_side_list, 0);
      requests.transfer(req2);
      requests.transfer(req3);
      if (comm_rank != 0) requests.wait();
      
      // filter sides here
      filter_sides<U>(
          ss_id,
          side_set_elem_list,
          side_set_side_list,
          cell_min,
          cell_max,
          element_global2local,
          side_id,
          side_cell,
          side_cell_face );
      
      start += chunk;
    } // while

    if (comm_rank == 0) {
      auto startt = MPI_Wtime();
      requests.wait();
      elapsed_time += MPI_Wtime() - startt;
      elapsed_counter++;
      std::cout << "Side set<" << ss_id << ">: Spent a total/average of " << elapsed_time 
        << " / " << elapsed_time / elapsed_counter << " secs waiting." << std::endl;
    }
  }

  return 0;
  
}

      
//============================================================================
//! \brief read the coordinates of the mesh from a file.
//! \param [in] exo_id  The exodus file id.
//! \return the vertex coordinates
//============================================================================
static int read_exo_side_sets(
    exo_reader_t & exo,
    const ex_init_params & exo_params,
    mpi_comm_t & comm,
    bool do_read,
    bool do_send,
    bool int64,
    size_t cell_min,
    size_t cell_max,
    const std::map<long_t, int_t> & cell_global2local,
    std::map<int_t, side_set_info_t> & side_sets,
    std::vector<int_t> & side_ids,
    std::vector<int_t> & side_cell,
    std::vector<int_t> & side_cell_face,
    std::vector<int_t> & side_offsets,
    std::vector<int_t> & sides)
{
  auto num_side_sets = exo_params.num_side_sets;
  if (!num_side_sets) return 0;

  // get the side set ids
  std::vector<int_t> ss_ids;

  if (do_read) {
    int ret;
    if (int64)
      ret = exo.read_side_set_ids<long long>(num_side_sets, ss_ids);
    else
      ret = exo.read_side_set_ids<int>(num_side_sets, ss_ids);
    if (ret) return ret;
  }
  if (do_send) broadcast_vector(comm, ss_ids, 0);

  // get the side set names
  std::vector<std::string> ss_names;
  if (do_read) {
    auto ret = exo.read_side_set_names(num_side_sets, ss_names);
    if (ret) return ret; 
  }
  else
    ss_names.resize(num_side_sets);

  // resize
  side_offsets.clear();
  side_offsets.reserve(num_side_sets+1);
  side_offsets.emplace_back(0);

  sides.clear();
  sides.reserve(num_side_sets);

  for (int i = 0; i < num_side_sets; i++){
    
    auto ss_id = ss_ids[i]-1;
    auto side_start = side_ids.size();
    
    // if no label, use the id
    if ( do_read && ss_names[i].empty() )
      ss_names[i] = std::to_string( ss_ids[i] ); 
    if (do_send) broadcast_string(comm, ss_names[i], 0);
  
    int ret;
    if (int64) {
      ret = read_exo_side_set<long long>(
          exo,
          comm,
          ss_ids[i],
          cell_min,
          cell_max,
          do_read,
          do_send,
          cell_global2local,
          side_ids,
          side_cell,
          side_cell_face);
    }
    else {
      ret = read_exo_side_set<int>(
          exo,
          comm,
          ss_ids[i],
          cell_min,
          cell_max,
          do_read,
          do_send,
          cell_global2local,
          side_ids,
          side_cell,
          side_cell_face);
    }
    if (ret) return ret;

    // keep track of the offsets
    auto side_end = side_ids.size();
    auto num_sides = side_end - side_start;
    if (num_sides) {
      side_offsets.emplace_back(side_end);
      sides.emplace_back(ss_id);
    }
    
    // if this side set is used on this rank
    side_sets.emplace(ss_id, side_set_info_t{ss_id, ss_names[i]});

  } // for side set 

  return 0;
}

//============================================================================
//! \brief read exodus faces
//============================================================================
template<typename U>
int read_exo_faces(
    exo_reader_t & exo,
    const ex_init_params & exo_params,
    mpi_comm_t & comm,
    bool do_read,
    bool do_send,
    const std::vector<ex_elem_blk_parm> & block_params,
    const crs_t<long_t, int_t> & cell_entities,
    int_t num_regions,
    const std::vector<int_t> & region_cell_offsets,
    const std::vector<shape_t> & cell_type,
    std::vector<int_t> & face_owner,
    crs_t<long_t, int_t> & face_vertices,
    crs_t<long_t, int_t> & cell_vertices,
    crs_t<long_t, int_t> & cell_faces)
{
  
  int ret;

  //----------------------------------------------------------------------------
  // read face ids
    
  auto num_face_blk = exo_params.num_face_blk;
  std::vector<int_t> face_blk_ids;
    
  if (do_read) {
    ret = exo.read_block_ids<U>(
        EX_FACE_BLOCK,
        num_face_blk,
        face_blk_ids);
    if (ret) {
      std::cout << "Error reading face block ids, received return status";
      std::cout << " of " << ret << std::endl;
      return ret;
    }
  }
  if (do_send) broadcast_vector(comm, face_blk_ids, 0);
  
  //----------------------------------------------------------------------------
  // read all face data
  
  std::map<long_t, int_t> face_global2local;
    
  std::map<std::vector<long_t>, long_t> sorted_face_vertices;

  size_t face_start = 0;
  
  for (int_t iblk=0; iblk<num_face_blk; ++iblk) {
    size_t num_faces;
    ret = read_exo_face_block<U>(
        exo,
        face_blk_ids[iblk],
        comm,
        do_read,
        do_send,
        cell_entities,
        num_regions,
        region_cell_offsets,
        cell_type,
        face_start,
        num_faces,
        face_global2local,
        face_vertices,
        sorted_face_vertices);
    if (ret) {
      std::cout << "Error reading face block, return status";
      std::cout << " of " << ret << std::endl;
      return ret;
    }
    face_start += num_faces;
  } // blocks

  // reset face owners
  face_owner.clear();
  face_owner.resize( face_vertices.size(), -1 );
  
  //----------------------------------------------------------------------------
  // Now build cell->faces and cell->vertices conn
  
  // first add all faces for fixed vertex entity types
  for (int_t r=0; r<num_regions; ++r) {
    auto region_start = region_cell_offsets[r];
    auto region_end = region_cell_offsets[r+1];
    auto elem_type = block_params[r].elem_type_val;

    //--- polyhedrons
    if (cell_type[region_start] == shape_t::polyhedron) {

      std::set<long_t> sorted_vs;
      std::vector<long_t> cell_vs, cell_fs;
      
      // loop over cells
      for (int_t c=region_start; c<region_end; ++c) {

        sorted_vs.clear();
        cell_vs.clear();
        cell_fs.clear();
        
        for (auto f : cell_entities.at(c)) {
          // map to local face id
          auto fid = face_global2local.at(f);
          // is this the first cell that uses this face
          if (face_owner[fid] == -1) face_owner[fid] = c;
          // add face to cell list
          cell_fs.emplace_back(fid);
          // find unique vertex list
          for (auto v : face_vertices.at(fid)) {
            auto ret = sorted_vs.emplace(v);
            if (ret.second) cell_vs.emplace_back(v);
          }
        }

        cell_vertices.push_back( cell_vs.begin(), cell_vs.end() );
        cell_faces.push_back( cell_fs.begin(), cell_fs.end() );

      } // cells

    }
    //--- fixed types
    else {

      // get side->vert info
      int table_len;
      int num_sides;
      int * side_vert_map;
      elem_side_table(elem_type, side_vert_map, table_len, num_sides);
      if (!table_len || !side_vert_map )
        THROW_ERROR("Face deduction not implemented for this block type");
            
      std::vector<long_t> cell_fs, face_vs, sorted_vs;
      cell_fs.reserve(table_len);
      face_vs.reserve(table_len);
      sorted_vs.reserve(table_len);

      // loop over cells
      for (int_t c=region_start; c<region_end; ++c) {

        cell_fs.clear();
        const auto & cell_vs = cell_entities.at(c);
        int_t num_verts = cell_vs.size();

        for (int s=0; s<num_sides; ++s) {
          
          face_vs.clear();
          for (int v=0; v<table_len; ++v) {
            auto id = side_vert_map[s*table_len + v];
            if (id < 1 || id > num_verts) continue;
            face_vs.emplace_back( cell_vs[id-1] );
          }

          sorted_vs.assign(face_vs.begin(), face_vs.end());
          std::sort( sorted_vs.begin(), sorted_vs.end() );
          auto ret = sorted_face_vertices.emplace(
              sorted_vs,
              sorted_face_vertices.size());
          if (ret.second) {
            face_vertices.push_back( face_vs.begin(), face_vs.end() );
            face_owner.emplace_back(c);
          }
          
          cell_fs.emplace_back(ret.first->second);
        } // faces

        cell_vertices.push_back( cell_vs.begin(), cell_vs.end() );
        cell_faces.push_back( cell_fs.begin(), cell_fs.end() );


      } // cells

    } //--- cell type
  } // regions

  return 0;
}

} //namespace

#endif

