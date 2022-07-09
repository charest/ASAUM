#include "partition.hpp"
#include "comm/comm_utils.hpp"
#include "comm/sub_comm.hpp"
#include "math/subdivide.hpp"
#include "utils/crs.hpp"
#include "utils/errors.hpp"


#ifdef HAVE_PARMETIS
#include <metis.h>
#include <parmetis.h>
#endif

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Create a distributed graph
////////////////////////////////////////////////////////////////////////////////
void graph(
    mpi_comm_t & comm,
    const std::vector<int_t> block_offsets,
    const std::vector< const crs_t<int_t, int_t> * > & cell2vert,
    const std::vector< const std::vector<long_t> * > & vert_local2global,
    int_t min_connections,
    std::vector<size_t> & distribution,
    crs_t<size_t, int_t> & graph)
{
  size_t rank = comm.rank();
  size_t size = comm.size();

  //const std::vector<int_t> block_counts,

  //----------------------------------------------------------------------------
  // Get cell partitioning
  //
  // Note: We create a global index space, but this may be different
  // from the global numbering of the original mesh.
  //----------------------------------------------------------------------------
  
  // determine how many blocks there are
  size_t num_blocks = cell2vert.size();
  size_t tot_blocks = block_offsets.back();

  // determine how many cells each block have
  size_t local_cells = 0;
  std::vector<size_t> block_sizes(num_blocks);
  for (size_t b=0; b<num_blocks; ++b) {
    block_sizes[b] = cell2vert[b]->size();
    local_cells += block_sizes[b];
  }

  // exchange the cell-block distribution 
  std::vector<int> recvdispls(block_offsets.begin(), block_offsets.end());
  std::vector<int> recvcounts(size);
  for(size_t r = 0; r < size; ++r)
    recvcounts[r] = recvdispls[r+1] - recvdispls[r];

  distribution.clear();
  distribution.resize(tot_blocks+1);
  auto dref = make_array_ref(&distribution[1], tot_blocks);
  comm.all_gatherv(
      block_sizes,
      dref,
      recvcounts,
      recvdispls);
  
  distribution[0] = 0;
  for (size_t b=0; b<tot_blocks; ++b)
    distribution[b+1] += distribution[b];

  // starting and ending cell indices
  auto block_start = block_offsets[rank];
  auto cells_start = distribution[block_start];

  //----------------------------------------------------------------------------
  // Create a naive initial distribution of the vertices
  //----------------------------------------------------------------------------

  // first figure out the maximum global id on this rank
  size_t max_global_vert_id{0};
  for (auto i : vert_local2global)
    for(auto v : *i)
      max_global_vert_id = std::max<size_t>(max_global_vert_id, v);

  // now the global max id
  size_t tot_verts{0};
  comm.all_reduce(max_global_vert_id, tot_verts, redop_t::max);
  tot_verts++;

  std::vector<size_t> vert_dist;
  subdivide(tot_verts, size, vert_dist);

  //----------------------------------------------------------------------------
  // Create the Local vertex-to-cell graph.
  //----------------------------------------------------------------------------

  // essentially vertex to cell connectivity
  std::map<size_t, std::vector<size_t>> vertex2cell;

  // Travel from the FROM_DIMENSION (cell) to the TO_DIMENSION (other)
  for (size_t b=0, cnt=0; b<num_blocks; ++b) {
    const auto & c2v = *cell2vert[b];
    const auto & v_l2g = *vert_local2global[b];
    size_t num_cells = c2v.size();
    const auto & offsets = c2v.offsets;
    const auto & indices = c2v.indices;
    for(size_t ic = 0; ic < num_cells; ++ic, ++cnt) {
      auto cell = cells_start + cnt;
      for(auto i = offsets[ic]; i < offsets[ic + 1]; ++i) {
        auto iv = indices[i];
        auto vertex = v_l2g[iv];
        vertex2cell[vertex].emplace_back(cell);
      }
    }
  }

  // remove duplicates
  auto remove_duplicates = [](auto & vertex2cell) {
    for(auto & vertex_cells : vertex2cell) {
      auto & cells = vertex_cells.second;
      std::sort(cells.begin(), cells.end());
      auto first = cells.begin();
      auto last = std::unique(first, cells.end());
      cells.erase(last, cells.end());
    }
  };

  remove_duplicates(vertex2cell);

  //----------------------------------------------------------------------------
  // Send vertex information to the deemed vertex-owner
  //----------------------------------------------------------------------------

  // We send the results for each vertex to their owner rank, which was
  // defined using the vertex subdivision above
  std::vector<int> sendcounts(size, 0);

  // count
  for(const auto & vs_pair : vertex2cell) {
    size_t global_id = vs_pair.first;
    auto r = owner(vert_dist, global_id);
    // we will be sending vertex id, number of cells, plus cell ids
    if(r != rank) {
      auto n = 2 + vs_pair.second.size();
      sendcounts[r] += n;
    }
  }

  // finish displacements
  std::vector<int> senddispls(size + 1);
  senddispls[0] = 0;
  for(size_t r = 0; r < size; ++r) {
    senddispls[r + 1] = senddispls[r] + sendcounts[r];
    sendcounts[r] = 0;
  }

  // fill buffers
  std::vector<size_t> sendbuf(senddispls[size]);

  for(const auto & vs_pair : vertex2cell) {
    size_t global_id = vs_pair.first;
    auto r = owner(vert_dist, global_id);
    // we will be sending vertex id, number of cells, plus cell ids
    if(r != rank) {
      // get offset
      auto offset = senddispls[r] + sendcounts[r];
      // populate data
      sendbuf[offset++] = global_id;
      sendbuf[offset++] = vs_pair.second.size();
      for(auto v : vs_pair.second)
        sendbuf[offset++] = v;
      // bump counters
      auto n = 2 + vs_pair.second.size();
      sendcounts[r] += n;
    }
  }

  comm.all_to_all(sendcounts, recvcounts);

  // how much info will we be receiving
  recvdispls[0] = 0;
  for(size_t r = 0; r < size; ++r)
    recvdispls[r + 1] = recvdispls[r] + recvcounts[r];
  std::vector<size_t> recvbuf(recvdispls[size]);

  // now send the actual vertex info
  comm.all_to_allv(sendbuf, sendcounts, senddispls, recvbuf, recvcounts,
    recvdispls);

  //----------------------------------------------------------------------------
  // Append received vertex information to the local vertex-to-cell graph
  //----------------------------------------------------------------------------

  // and add the results to our local list
  std::map<size_t, std::vector<size_t>> vertex2rank;

  for(size_t r = 0; r < size; ++r) {
    for(auto i = recvdispls[r]; i < recvdispls[r + 1];) {
      // get vertex
      auto vertex = recvbuf[i];
      ++i;
      // keep track of ranks that share this vertex
      vertex2rank[vertex].emplace_back(r);
      assert(i < recvdispls[r + 1]);
      // unpack cell neighbors
      auto n = recvbuf[i];
      ++i;
      for(size_t j = 0; j < n; ++j) {
        assert(i < recvdispls[r + 1]);
        auto cell = recvbuf[i];
        ++i;
        // might not already be there
        vertex2cell[vertex].emplace_back(cell);
      }
    }
  }

  // remove duplicates
  remove_duplicates(vertex2cell);

  //----------------------------------------------------------------------------
  // Send back the final results for shared vertices
  //----------------------------------------------------------------------------

  // now perpare to send results back

  // count send buffer size
  std::fill(sendcounts.begin(), sendcounts.end(), 0);
  for(const auto & vertex_pair : vertex2rank) {
    auto vertex = vertex_pair.first;
    for(auto r : vertex_pair.second) {
      if(r != rank) { // should always enter anyway!
        // better already be there
        const auto & cells = vertex2cell.at(vertex);
        // we will be sending vertex id, number of cells, plus cell ids
        auto n = 2 + cells.size();
        sendcounts[r] += n;
      }
    }
  }

  // finish displacements
  senddispls[0] = 0;
  for(size_t r = 0; r < size; ++r)
    senddispls[r + 1] = senddispls[r] + sendcounts[r];

  // resize send buffer
  sendbuf.clear();
  sendbuf.resize(senddispls[size]);

  // now fill buffer
  std::fill(sendcounts.begin(), sendcounts.end(), 0);
  for(const auto & vertex_pair : vertex2rank) {
    auto vertex = vertex_pair.first;
    for(auto r : vertex_pair.second) {
      if(r != rank) { // should always enter anyway!
        auto j = senddispls[r] + sendcounts[r];
        // better already be there
        const auto & cells = vertex2cell.at(vertex);
        // we will be sending vertex id, number of cells, plus cell ids
        sendbuf[j++] = vertex;
        sendbuf[j++] = cells.size();
        for(auto c : cells)
          sendbuf[j++] = c;
        // increment send counter
        sendcounts[r] = j - senddispls[r];
      }
    }
  }

  // send counts
  std::fill(recvcounts.begin(), recvcounts.end(), 0);
  comm.all_to_all(sendcounts, recvcounts);

  // how much info will we be receiving
  recvdispls[0] = 0;
  for(size_t r = 0; r < size; ++r)
    recvdispls[r + 1] = recvdispls[r] + recvcounts[r];

  // how much info will we be receiving
  recvbuf.clear();
  recvbuf.resize(recvdispls.at(size));

  // now send the final vertex info back
  comm.all_to_allv(sendbuf, sendcounts, senddispls, recvbuf, recvcounts,
    recvdispls);

  //----------------------------------------------------------------------------
  // Append received vertex information to the local vertex-to-cell graph
  //----------------------------------------------------------------------------

  for(size_t r = 0; r < size; ++r) {
    for(auto i = recvdispls[r]; i < recvdispls[r + 1];) {
      // get vertex
      auto vertex = recvbuf[i];
      ++i;
      // keep track of ranks that share this vertex
      vertex2rank[vertex].emplace_back(r);
      assert(i < recvdispls[r + 1]);
      // unpack cell neighbors
      auto n = recvbuf[i];
      ++i;
      for(size_t j = 0; j < n; ++j) {
        assert(i < recvdispls[r + 1]);
        auto cell = recvbuf[i];
        ++i;
        // might not already be there
        vertex2cell[vertex].emplace_back(cell);
      }
    }
  }

  // remove duplicates
  remove_duplicates(vertex2cell);

  //----------------------------------------------------------------------------
  // Invert for final cell-to-cell graph
  //----------------------------------------------------------------------------
  
  graph.clear();
  auto & graph_offsets = graph.offsets;
  auto & graph_indices = graph.indices;

  // Set the first offset (always zero).
  graph_offsets.reserve(local_cells);
  //graph.indices.reserve(n * num_cells); // guess
  graph_offsets.push_back(0);

  std::map<size_t, size_t> cell_counts;

  for (size_t b=0, cnt=0; b<num_blocks; ++b) {
    
    const auto & c2v = *cell2vert[b];
    const auto & cell_offsets = c2v.offsets;
    const auto & cell_indices = c2v.indices;
    
    size_t num_cells = c2v.size();
    
    const auto & v_l2g = *vert_local2global[b];

    for(size_t ic = 0; ic < num_cells; ++ic, ++cnt) {
      auto cell = cells_start + cnt;

      cell_counts.clear();

      // iterate over vertices
      auto start = cell_offsets[ic];
      auto end = cell_offsets[ic + 1];
      for(auto i = start; i < end; ++i) {
        auto iv = cell_indices[i];
        auto vertex = v_l2g[iv];
        // now count attached cells
        for(auto other : vertex2cell.at(vertex)) {
          // new entries are initialized to zero
          if(other != cell)
            cell_counts[other] += 1;
        }
      }

      // now add results
      int num_connections{0};

      for(auto count : cell_counts) {
        if(count.second >= static_cast<size_t>(min_connections))
        {
          graph_indices.emplace_back(count.first);
          num_connections++;
        }
      }

      graph_offsets.emplace_back(graph_offsets.back() + num_connections);

    } // cells
  } // blocks

} // make_dcrs

////////////////////////////////////////////////////////////////////////////////
/// Naiive partitioner
////////////////////////////////////////////////////////////////////////////////
void naive(
    mpi_comm_t & comm,
    const std::vector<int_t> & block_offsets,
    const std::vector<size_t> & distribution,
    int_t tot_parts,
    std::vector<int_t> & parts)
{
  size_t tot_ents = distribution.back();
  std::vector<size_t> newdist;
  subdivide(tot_ents, tot_parts, newdist);

  auto rank = comm.rank();

  auto block_start = block_offsets[rank];
  auto block_end = block_offsets[rank+1];
  size_t num_ents = distribution[block_end] - distribution[block_start];
  
  parts.clear();
  parts.reserve(num_ents);

  for (auto b=block_start; b<block_end; ++b) {
    size_t cells_start = distribution[b];
    size_t num_ents = distribution[b+1] - cells_start;
    for (size_t i=0, r=0; i<num_ents; ++i) {
      auto global_id = cells_start + i;
      while (global_id>=newdist[r+1]) {r++;}
      parts.emplace_back(r);
    }
  }

}


////////////////////////////////////////////////////////////////////////////////
/// METIS partitioner
////////////////////////////////////////////////////////////////////////////////
void metis(
    mpi_comm_t & comm,
    const std::vector<int_t> & block_offsets,
    const std::vector<size_t> & distribution,
    const crs_t<size_t, int_t> & graph,
    int_t tot_parts,
    std::vector<int_t> & parts_out)
{
#ifdef HAVE_PARMETIS
  parts_out.clear();

  // only populated partitions need participate
  size_t num_ents = graph.size();
  sub_comm_t subcomm(comm, num_ents);
  if (!num_ents) return;

  // comm stats
  int_t comm_size  = subcomm.comm().size();

  // subdivide the partition distributions evenly
  std::vector<int_t> part_dist;
  subdivide(tot_parts, comm_size, part_dist);
  
  // how many pieces is this rank to create
  size_t comm_rank = subcomm.comm().rank();
  auto num_pieces = part_dist[comm_rank+1] - part_dist[comm_rank];
  
  // Had a piece, now there is none, this woud require redistribution, and is not allowed 
  if (num_pieces == 0) {
    THROW_ERROR("Local subdivision not possible for this configuration."
        << " Either change how the original partitions are distributed,"
        << " or use a parallel partitioner");
  }
  // no partitioning necessary
  else if (num_pieces == 1 ) {
    parts_out.resize(num_ents);
    std::fill(parts_out.begin(), parts_out.end(), part_dist[comm_rank]);
    return;
  }
  
  
  // Get the dCRS information using ParMETIS types.  Also convert the global graph
  // to a local one
  size_t world_rank = comm.rank();
  auto blkstart = block_offsets[world_rank];
  auto blkend = block_offsets[world_rank+1];
  
  auto vtxstart = distribution[blkstart];
  auto vtxend = distribution[blkend];

  const auto & graph_offsets = graph.offsets;
  const auto & graph_indices = graph.indices;

  std::vector<idx_t> xadj, adjncy;
  xadj.reserve(graph_offsets.size());
  adjncy.reserve(graph_indices.size());
      
  xadj.emplace_back(0);
  for (size_t i=0; i<num_ents; ++i) {
    for (auto j=graph_offsets[i]; j<graph_offsets[i+1]; ++j) {
      auto id = graph_indices[j];
      if (id >= vtxstart && id < vtxend)
        adjncy.emplace_back(id - vtxstart);
    }
    xadj.emplace_back(adjncy.size());
  }

  // set inputs
  idx_t ncon = 1;
  std::vector<::real_t> tpwgts(ncon * num_pieces, 1.0 / num_pieces);
  
  // We may need to expose some of the ParMETIS configuration options.
  std::vector<::real_t> ubvec(ncon, 1.05);
  idx_t edgecut;

  std::vector<idx_t> part(num_ents);
  
  // Actual call to ParMETIS.
  idx_t nvtx = num_ents; 
  idx_t npart = num_pieces;
  int result = METIS_PartGraphKway(&nvtx, &ncon, &xadj[0], &adjncy[0],
    nullptr, nullptr, nullptr, &npart, &tpwgts[0],
    ubvec.data(), nullptr, &edgecut, &part[0]);
  if(result != METIS_OK)
    THROW_ERROR("Metis' METIS_PartGraphKway failed!");
  
  parts_out.assign(part.begin(), part.end());
  for (auto & p : parts_out) p += part_dist[comm_rank];

#else

  THROW_ERROR("Rebuild with METIS to enable METIS partitioning.");

#endif
}

////////////////////////////////////////////////////////////////////////////////
/// PARMETIS graph partitioner
////////////////////////////////////////////////////////////////////////////////
void parmetis_graph(
    mpi_comm_t & comm,
    const std::vector<int_t> & block_offsets,
    const std::vector<size_t> & distribution,
    const crs_t<size_t, int_t> & graph,
    int_t tot_parts,
    std::vector<int_t> & parts_out,
    const std::vector<real_t> & xyz = {})
{
#ifdef HAVE_PARMETIS
  parts_out.clear();

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  auto num_ents = graph.size();
  sub_comm_t subcomm(comm, num_ents);
  if (num_ents == 0) return;

  // Call ParMETIS partitioner.
  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon = 1;
  std::vector<::real_t> tpwgts(ncon * tot_parts, 1.0 / tot_parts);
  
  // We may need to expose some of the ParMETIS configuration options.
  std::vector<::real_t> ubvec(ncon, 1.05);
  idx_t options[3] = {0, 0, 0};
  idx_t edgecut;
  std::vector<idx_t> part(graph.size());
  
  // Get the dCRS information using ParMETIS types.
  std::vector<idx_t> xadj(graph.offsets.begin(), graph.offsets.end());
  std::vector<idx_t> adjncy(graph.indices.begin(), graph.indices.end());
  
  // colapse distribution
  size_t world_size = comm.size();
  size_t comm_size = subcomm.comm().size();

  std::vector<idx_t> vtxdist(comm_size+1);
  vtxdist[0] = 0;

  for (size_t r=0, i=0; r<world_size; ++r) {
    auto blkstart = block_offsets[r];
    auto blkend = block_offsets[r+1];
    auto n = distribution[blkend] - distribution[blkstart];
    if (n>0) {
      vtxdist[i+1] = vtxdist[i] + n;
      i++;
    }
  }
  
  // Actual call to ParMETIS.
  idx_t npart = tot_parts;
  MPI_Comm mpi_comm = subcomm.comm().comm();

  if (xyz.empty()) {
    int result = ParMETIS_V3_PartKway(&vtxdist[0], &xadj[0], &adjncy[0],
        nullptr, nullptr, &wgtflag, &numflag, &ncon, &npart, &tpwgts[0],
        ubvec.data(), options, &edgecut, &part[0], &mpi_comm);
    if(result != METIS_OK)
      THROW_ERROR("Parmetis' ParMETIS_V3_PartKway failed!");
  }
  else {
    std::vector<::real_t> _xyz(xyz.begin(), xyz.end());
    idx_t ndims = xyz.size() / num_ents;
    int result = ParMETIS_V3_PartGeomKway(&vtxdist[0], &xadj[0], &adjncy[0],
        nullptr, nullptr, &wgtflag, &numflag, &ndims, _xyz.data(), &ncon,
        &npart, &tpwgts[0], ubvec.data(), options, &edgecut, &part[0],
        &mpi_comm);
    if(result != METIS_OK)
      THROW_ERROR("Parmetis' ParMETIS_V3_PartGeomKway failed!");
  }
  
  parts_out.assign(part.begin(), part.end());

#else

  THROW_ERROR("Rebuild with ParMETIS to enable ParMETIS partitioning.");

#endif
}


////////////////////////////////////////////////////////////////////////////////
/// PARMETIS geom partitioner
////////////////////////////////////////////////////////////////////////////////
void parmetis_geom(
    mpi_comm_t & comm,
    const std::vector<int_t> & block_offsets,
    const std::vector<size_t> & distribution,
    const std::vector<real_t> & xyz,
    int_t tot_parts,
    std::vector<int_t> & parts_out)
{
#ifdef HAVE_PARMETIS
  parts_out.clear();

  auto world_rank = comm.rank();

  auto block_start = block_offsets[world_rank];
  auto block_end = block_offsets[world_rank+1];
  size_t num_ents = distribution[block_end] - distribution[block_start];
  
  std::vector<idx_t> vtxdist(distribution.begin(), distribution.end());

  // create a subcommunicator
  sub_comm_t subcomm(comm, num_ents);
  if (num_ents == 0) return;
  
  int_t comm_size = subcomm.comm().size();
  if (tot_parts != comm_size) {
    THROW_ERROR("Unfortunately, ParMETIS_V3_PartGeom requires nparts == comm_size."
        << " " << comm_size << " ranks are tryinig to compute a partition of "
        << "size " << tot_parts << ".");
  }
  
  std::vector<idx_t> part;
  part.reserve(num_ents);
  for (auto b=block_start; b<block_end; ++b) {
    size_t n = distribution[b+1] - distribution[b];
    part.insert(part.end(), n, b);
  }

  // Actual call to ParMETIS.
  idx_t ndims = xyz.size() / num_ents;
  std::vector<::real_t> _xyz(xyz.begin(), xyz.end());
  MPI_Comm mpi_comm = subcomm.comm().comm();
  int result = ParMETIS_V3_PartGeom(&vtxdist[0], &ndims, _xyz.data(), part.data(), &mpi_comm);
  if(result != METIS_OK)
    THROW_ERROR("Parmetis' ParMETIS_V3_PartGeom failed!");
  
  parts_out.assign(part.begin(), part.end());

#else

  THROW_ERROR("Rebuild with ParMETIS to enable ParMETIS partitioning.");

#endif
}


////////////////////////////////////////////////////////////////////////////////
/// Partition a mesh
////////////////////////////////////////////////////////////////////////////////
void partition(
    mpi_comm_t & comm,
    const std::vector<int_t> & block_offsets,
    const std::vector<size_t> & distribution,
    const crs_t<size_t, int_t> & grph,
    part_alg_t part_alg,
    const std::vector<real_t> & cell_xyz,
    int_t tot_parts,
    std::vector<int_t> & parts)
{
  auto is_root = comm.is_root();

  // now partition
  if (part_alg == part_alg_t::parmetis_kway) {
    if (is_root) std::cout << "msh> Partitioning via ParMETIS KWAY." << std::endl;
    parmetis_graph(
        comm,
        block_offsets,
        distribution,
        grph,
        tot_parts,
        parts);
  }
  else if (part_alg == part_alg_t::parmetis_geomkway) {
    if (is_root) std::cout << "msh> Partitioning via ParMETIS Geom-KWAY." << std::endl;
    if (cell_xyz.empty())
      THROW_ERROR("You need to provide coordinates to use"
          << " geometry-based partitioners.");
    parmetis_graph(
        comm,
        block_offsets,
        distribution,
        grph,
        tot_parts,
        parts,
        cell_xyz);
  }
  else if (part_alg == part_alg_t::parmetis_geom) {
    if (is_root) std::cout << "msh> Partitioning via ParMETIS Geom." << std::endl;
    if (cell_xyz.empty())
      THROW_ERROR("You need to provide coordinates to use"
          << " geometry-based partitioners.");
    parmetis_geom(
        comm,
        block_offsets,
        distribution,
        cell_xyz,
        tot_parts,
        parts);
  }
  else if (part_alg == part_alg_t::metis) {
    if (is_root) std::cout << "msh> Partitioning via METIS." << std::endl;
    metis(
        comm,
        block_offsets,
        distribution,
        grph,
        tot_parts,
        parts);
  }
  else if (part_alg == part_alg_t::naive) {
    if (is_root) std::cout << "msh> Partitioning via naive subdivision." << std::endl;
    naive(
        comm,
        block_offsets,
        distribution,
        tot_parts,
        parts);
  }
  else {
    THROW_ERROR("Unknown partitioner selected.");
  }
}

} // namespace
