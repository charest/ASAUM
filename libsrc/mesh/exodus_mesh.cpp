#include "exodus_mesh.hpp"
#include "exodus_mesh_utils.hpp"
#include "mesh_unstruct_block.hpp"
#include "comm/comm_map.hpp"
#include "comm/comm_utils.hpp"
#include "comm/entity_info.hpp"
#include "comm/sub_comm.hpp"
#include "math/subdivide.hpp"
#include "lua/lua_utils.hpp"
#include "utils/cast.hpp"

#include <unordered_map>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Box mesh constructor
////////////////////////////////////////////////////////////////////////////////
exodus_mesh_t::exodus_mesh_t(
    lua_result_t input,
    mpi_comm_t & comm,
    const std::vector<int_t> & parts,
    const std::string & filename,
    distribution_t distribute,
    part_alg_t part_alg,
    bool repart) :
  mesh_t(input, comm, false, parts),
  filename_(filename),
  distribution_(distribute),
  part_alg_(part_alg),
  repart_(repart)
{
  // get communicator info
  auto is_root = comm.is_root();

  if (parts.size() > 1) {
    if (is_root) {
      std::cout << "Unstructured meshes can only be partitioned by in ";
      std::cout << "one dimension.  You provided " << parts.size();
      std::cout << " values.  Only one allowed.";
    }
    comm.exit(consts::failure);
  }

  read_and_partition();

  number();

}

////////////////////////////////////////////////////////////////////////////////
/// Load the mesh
////////////////////////////////////////////////////////////////////////////////
void exodus_mesh_t::read_and_partition()
{
  // get communicator info
  auto is_root = comm_.is_root();
  auto comm_size = comm_.size();
  auto comm_rank = comm_.rank();

  if (is_root)
    std::cout << "msh> Loading exodus mesh: " << filename_ << std::endl;

  //----------------------------------------------------------------------------
  // Initial setup

  // Figure out how many initial partitions we will have
  bool serial = true;
  int_t tot_parts = comm_size;
  std::vector<int_t> part_ids = {comm_rank};
  std::vector<std::string> filenames = {filename_};
    
  // figure how many initial partitions we will have
  auto ext = file_extension(filename_);

  //--- Multiple files
  if (is_signed_integer(ext)) {
    serial = false;
    tot_parts = std::stoi(ext);
  
    if (tot_parts<1) {
      std::cout << "Error processing group of exodus files '" << filename_;
      std::cout << "', looks like there are no files in the group: ";
      std::cout << "'" << ext << "'" << std::endl;
      comm_.exit(consts::failure);
    }
  
    // distribute the initial partitions among processors
    if (distribution_ == distribution_t::sequential) {
      if (is_root) std::cout << "msh> Sequential file distribution." << std::endl;
      part_ids = distribute_sequentially(tot_parts, comm_size, comm_rank);
    }
    else if (distribution_ == distribution_t::cyclic) {
      if (is_root) std::cout << "msh> Cyclic file distribution." << std::endl;
      part_ids = distribute_cyclic(tot_parts, comm_size, comm_rank);
    }
    else if (distribution_ == distribution_t::hostname) {
      if (is_root) std::cout << "msh> Hostname file distribution." << std::endl;
      part_ids = distribute_key(tot_parts, comm_, comm_.hostname());
    }
 
    // how many digits in padded number
    auto ndigits = ext.size();
    // figure out this processors filenames
    auto num_parts = part_ids.size();
    filenames.resize(num_parts);
    for (size_t f=0; f<num_parts; ++f)
      filenames[f] = filename_ + "." + zero_padded(part_ids[f], ndigits);
  }
  
  // read the mesh
  read(comm_, serial, filenames, part_ids);
  if (is_root) std::cout << "msh> Finished reading." << std::endl;

  // some ranks might not know haw many dimensions there are, reduce them
  int_t max_dims;
  comm_.all_reduce(num_dims_, max_dims, redop_t::max);
  
  if (num_dims_ != 0 && num_dims_ != max_dims) {
    std::cout << "Some ranks have differing number of dimensions.  Check";
    std::cout << " file consistency." << std::endl;
    comm_.exit(consts::failure); 
  }

  num_dims_ = max_dims;

  // build the offsets
  int_t num_parts = num_blocks();
  global_offsets(comm_, num_parts, block_offsets_); 

  auto tot_blks = block_offsets_.back();
  if (!partitioning_.empty()) tot_parts = partitioning_[0];

  //----------------------------------------------------------------------------
  // partition the mesh
  if ((serial && tot_blks>1) || (tot_blks != tot_parts) || repart_)
  {

    if (is_root) std::cout << "msh> Needs partitioning." << std::endl;

    // get the cell connectivity
    auto tot_cells = num_cells();

    // if this needs cell coordinnates
    std::vector<real_t> cell_midpoints;
    
    if (part_alg_ == part_alg_t::parmetis_geom ||
        part_alg_ == part_alg_t::parmetis_geomkway)
    {
      cell_midpoints.reserve(num_dims_ * tot_cells);
      for (int_t b=0; b<num_parts; ++b) {
        auto this_block = dynamic_cast<mesh_unstruct_block_t*>(blocks_[b]);
        this_block->cell_midpoints(cell_midpoints);
      }
    }

    std::vector<size_t> dist;
    crs_t<size_t, int_t> grph;
    graph(num_dims_, dist, grph);

    std::vector<int_t> partitioning;
    partition(
        comm_,
        block_offsets_,
        dist,
        grph,
        part_alg_,
        cell_midpoints,
        tot_parts,
        partitioning);
    if (is_root) std::cout << "msh> Finished partitioning." << std::endl;

    // now migrate the entities to their respective ranks
    if (is_root) std::cout << "msh> Migrating the mesh." << std::endl;
    std::vector<int_t> part_dist;
    subdivide(tot_parts, comm_size, part_dist);
    migrate(partitioning, part_dist, part_ids); 
    if (is_root) std::cout << "msh> Finished migrating." << std::endl;

  } // partition
  
  //----------------------------------------------------------------------------
  // Post migration setup
  
  // create the roots
  tot_blks = tot_blocks();
  reserve_roots(tot_blks);
  for (int_t i=0; i<tot_blks; ++i)
    add_root( std::make_unique<mesh_unstruct_node_t>(num_dims_, i) );
  
  auto block_start = block_offsets_[comm_rank];
  auto block_end = block_offsets_[comm_rank+1];
  for (int_t b=0, bid=block_start; bid<block_end; ++b, ++bid) {
    auto this_block = dynamic_cast<mesh_unstruct_block_t*>(blocks_[b]);
    root(bid)->set_block(this_block);
    this_block->set_num_vertices();
    this_block->set_num_cells();
    if (num_dims_>1) this_block->set_num_entities(num_dims_-1);
    this_block->install_boundaries(boundaries_);
  }

}

////////////////////////////////////////////////////////////////////////////////
/// add a side set
////////////////////////////////////////////////////////////////////////////////
void exodus_mesh_t::add_side_set(const side_set_info_t & side_info)
{
  auto it = side_sets_.emplace(side_sets_.end(), side_info);
  side_id_map_[side_info.id] = it;
  side_name_map_[side_info.label] = it;
}

////////////////////////////////////////////////////////////////////////////////
/// Box mesh constructor
////////////////////////////////////////////////////////////////////////////////
void exodus_mesh_t::read(
    mpi_comm_t & comm,
    bool serial,
    const std::vector<std::string> & filenames,
    const std::vector<int_t> & part_ids)
{
 
  // create the communicator
  auto subcomm = std::make_unique<mpi_comm_t>(comm, part_ids.size());
  
  // if there is no files to read, return
  if (part_ids.empty()) return;

  // get communicator info
  auto is_root = subcomm->is_root();
  auto comm_size = subcomm->size();
  auto comm_rank = subcomm->rank();
    
  bool do_read = (serial && is_root) ||  (!serial);
  bool do_send = serial && (comm_size>1);

  exo_reader_t exo;

  int_t num_parts = part_ids.size();
  blocks_.reserve(num_parts);
    
  crs_t<long_t, int_t> cell_entities;

  size_t tot_cells = 0;

  std::map<int_t, side_set_info_t> side_sets;
  
  for (int_t p=0; p<num_parts; ++p) {

    //--------------------------------------------------------------------------
    // open the exodus file
    auto ret = exo.open(filenames[p]);
    if (ret) {
      std::cout << "Error opening file '" << filenames[p];
      std::cout << "', received return status";
      std::cout << " of " << ret << std::endl;
      subcomm->exit(consts::failure);
    }
    
    // get the initialization parameters
    ex_init_params exo_params;
    if (do_read) {
      ret = exo.read_params(exo_params);
      if (ret) {
        std::cout << "Error getting exodus params, received return status";
        std::cout << " of " << ret << std::endl;
        subcomm->exit(consts::failure);
      }
    }
    if (do_send) subcomm->broadcast(exo_params, 0);
    
    // create a new block
    num_dims_ = exo_params.num_dim;
    auto part = new mesh_unstruct_block_t(num_dims_);

    // check the integer type used in the exodus file
    bool int64;
    if (do_read) int64 = exo.is_int64();
    if (do_send) subcomm->broadcast(int64, 0);
    
    //--------------------------------------------------------------------------
    // Figure out partitioning 
    
    // figure out the number of cells for this rank.
    // if this is reading a bunch of exodus files, there is
    // nothing to do
    size_t num_cells = exo_params.num_elem;
    size_t cell_min{0}, cell_max{num_cells-1};
    std::vector<size_t> cell_partitioning = {cell_min, cell_max+1};

    if (serial) {
      cell_partitioning.clear();
      subdivide( num_cells, comm_size, cell_partitioning );

      cell_min = cell_partitioning[comm_rank];
      cell_max = cell_partitioning[comm_rank+1] - 1;
      num_cells = cell_max - cell_min + 1;
    }

    // keep track of how many cells are read
    tot_cells += num_cells;

    //--------------------------------------------------------------------------
    // element blocks
    
    // offsets always have initial zero
    cell_entities.clear();
    cell_entities.offsets.emplace_back( 0 );
    size_t cell_counter{0};

    // the number of blocks and some storage for block ids
    auto num_elem_blk = exo_params.num_elem_blk;
    std::vector<int_t> elem_blk_ids;

    // get the element block ids
    if (do_read) {
      int ret;
      if (int64)
        ret = exo.read_block_ids<long long>(EX_ELEM_BLOCK, num_elem_blk,
            elem_blk_ids);
      else
        ret = exo.read_block_ids<int>(EX_ELEM_BLOCK, num_elem_blk,
            elem_blk_ids);
      if (ret) {
        std::cout << "Error reading block ids, received return status";
        std::cout << " of " << ret << std::endl;
        subcomm->exit(consts::failure);
      }
    }
    if (do_send) broadcast_vector(*subcomm, elem_blk_ids, 0);
    
    // read block info
    std::vector<ex_elem_blk_parm> elem_blk_params(elem_blk_ids.size());
    if (do_read) {
      auto ret = exo.read_block_params(EX_ELEM_BLOCK, elem_blk_ids, elem_blk_params);
      if (ret) {
        std::cout << "Error reading block params, received return status";
        std::cout << " of " << ret << std::endl;
        subcomm->exit(consts::failure);
      }
    }
    if (do_send) comm_.broadcast(elem_blk_params, 0);

    // read each block
    for (int iblk = 0; iblk < num_elem_blk; iblk++) {

      shape_t block_type;
      auto cell_start = cell_entities.size();
      auto blk_id = elem_blk_ids[iblk];

      // 64  bit version
      if (int64) 
        ret = read_exo_block<long long>(
            exo,
            blk_id,
            cell_min,
            cell_max,
            cell_partitioning,
            *subcomm,
            do_read,
            do_send,
            cell_counter,
            cell_entities,
            block_type);
      else 
        ret = read_exo_block<int>(
            exo,
            blk_id,
            cell_min,
            cell_max,
            cell_partitioning,
            *subcomm,
            do_read,
            do_send,
            cell_counter,
            cell_entities,
            block_type);
      if (ret) {
        std::cout << "Error reading blocks, received return status";
        std::cout << " of " << ret << std::endl;
        subcomm->exit(consts::failure);
      }
      
      if (block_type == shape_t::empty) continue;

      // add the block type
      size_t num_added = cell_entities.size() - cell_start;
      part->cell_type_.insert(part->cell_type_.end(), num_added, block_type);

      // add the regions
      auto region_id = blk_id-1;
      part->cell_region_.insert(part->cell_region_.end(), num_added, region_id);
      part->add_region(region_id, num_added);

    } // blocks

    // check some assertions
    if (ret) {
      std::cout << "Mismatch in read blocks. Expected " << num_cells;
      std::cout << " got " << cell_entities.size() << std::endl;
      subcomm->exit(consts::failure);
    }
    
    // create the local to global cell mapping
    auto & cell_local2global_ = part->local2global_[num_dims_];
    auto & cell_global2local_ = part->global2local_[num_dims_];
    
    ret = read_exo_cell_ids(
        exo,
        serial,
        int64,
        cell_min,
        num_cells,
        cell_local2global_,
        cell_global2local_);
    if (ret) {
      std::cout << "Error reading global ids, return status";
      std::cout << " of " << ret << std::endl;
      subcomm->exit(consts::failure);
    }

    //--------------------------------------------------------------------------
    // Read faces
    
    crs_t<long_t, int_t> cell_vertices, cell_faces;
    crs_t<long_t, int_t> face_vertices;

    if (num_dims_ < 3) {

      std::swap(cell_vertices, cell_entities);

    }
    else {
    
      if (int64)
        ret = read_exo_faces<long long>(
            exo,
            exo_params,
            *subcomm,
            do_read,
            do_send,
            elem_blk_params,
            cell_entities,
            part->regions_.size(),
            part->region_cell_offsets_,
            part->cell_type_,
            part->face_owner_,
            face_vertices,
            cell_vertices,
            cell_faces);
      else
        ret = read_exo_faces<int>(
            exo,
            exo_params,
            *subcomm,
            do_read,
            do_send,
            elem_blk_params,
            cell_entities,
            part->regions_.size(),
            part->region_cell_offsets_,
            part->cell_type_,
            part->face_owner_,
            face_vertices,
            cell_vertices,
            cell_faces);
      if (ret) {
        std::cout << "Error reading faces, return status";
        std::cout << " of " << ret << std::endl;
        subcomm->exit(consts::failure);
      }

    } // num dims
        
    //--------------------------------------------------------------------------
    // Read vertices

    auto & vertex_local2global_ = part->local2global_[0];
    auto & vertex_global2local_ = part->global2local_[0];
    
    // get the vertex maps
    if ( serial ) {

      ret = read_exo_point_coords_serial(
          exo,
          exo_params,
          *subcomm,
          do_send,
          cell_vertices,
          vertex_global2local_,
          vertex_local2global_,
          part->vertex_coords_);

      // convert element conectivity to local ids 
      for (auto & v : cell_vertices.indices)
        v = vertex_global2local_.at(v);
      for (auto & v : face_vertices.indices)
        v = vertex_global2local_.at(v);

    } // serial
    // parallel
    else {

      ret = read_exo_point_coords_parallel(
          exo,
          exo_params,
          int64,
          vertex_global2local_,
          vertex_local2global_,
          part->vertex_coords_);
      
    } // parallel
    
    if (ret) {
      std::cout << "Error reading vertex coords, return status";
      std::cout << " of " << ret << std::endl;
      subcomm->exit(consts::failure);
    }
      
    //--------------------------------------------------------------------------
    // Copy over connectivity
    
    // get reference to cell->vertex connectivity
    auto & cell_vertices_ = part->connectivity_[{num_dims_, 0}];
    // copy over
    cell_vertices_.assign( cell_vertices );

    if (!face_vertices.empty())
      part->connectivity_[{num_dims_-1, 0}].assign(face_vertices);

    if (!cell_faces.empty())
      part->connectivity_[{num_dims_, num_dims_-1}].assign(cell_faces);

    //--------------------------------------------------------------------------
    // read side sets
    ret = read_exo_side_sets(
        exo,
        exo_params,
        *subcomm,
        do_read,
        do_send,
        int64,
        cell_min,
        cell_max,
        cell_global2local_,
        side_sets,
        part->side_id_,
        part->side_cell_,
        part->side_cell_face_,
        part->side_set_offsets_,
        part->side_sets_);
    if (ret) {
      std::cout << "Error reading side sets, return status";
      std::cout << " of " << ret << std::endl;
      subcomm->exit(consts::failure);
    }

    //--------------------------------------------------------------------------
    // close the file
    if (do_read) exo.close();

    // any post setup
    if (num_dims_>1) part->sort_face_vertices();
    part->build_cell_sides();
    part->set_num_vertices();
    part->set_num_cells();
    if (num_dims_>1) part->set_num_entities(num_dims_-1);
    
    // move the partition
    blocks_.emplace_back( part );

  } // files
     
  //----------------------------------------------------------------------------
  // Finish side sets

  // aglomerate side set info
  side_set_info_t::broadcast( side_sets, *subcomm );
  for (const auto & si : side_sets)
    add_side_set(si.second);

}

////////////////////////////////////////////////////////////////////////////////
/// migrate a mesh
////////////////////////////////////////////////////////////////////////////////
void exodus_mesh_t::migrate(
    const std::vector<int_t> & cell_parts,
    const std::vector<int_t> & part_dist,
    const std::vector<int_t> & part_ids)
{

  size_t world_size = comm_.size();
  size_t world_rank = comm_.rank();
  
  //----------------------------------------------------------------------------
  // Pack information to migrate
  //----------------------------------------------------------------------------
  
  //------------------------------------
  // create storage for counters

  auto num_blks = num_blocks();

  auto tot_parts = part_dist.back();
  auto num_parts = part_dist[world_rank+1] - part_dist[world_rank];

  std::vector<byte_t> tempbuf;
  std::vector<int> sendcounts(tot_parts, 0);
  std::vector<int> localcounts(num_parts, 0);
  std::vector<int_t> block_erase_counts(num_blks, 0);

  // mapping between global and local partition ids
  std::map<int_t, int_t> part_map;
  for (auto b=part_dist[world_rank]; b<part_dist[world_rank+1]; ++b)
    part_map[b] = b - part_dist[world_rank];
        
  std::vector<int_t> rank_owner;
  invert_owner(part_dist, rank_owner);

  //------------------------------------
  // Count sizes
  for (int_t b=0, cnt=0; b<num_blks; ++b) {

    auto this_block = dynamic_cast<mesh_unstruct_block_t*>(blocks_[b]);
    size_t num_elements = this_block->local2global_.at(num_dims_).size();

    auto & erase_counter = block_erase_counts[b];

    for(size_t local_id(0); local_id < num_elements; ++local_id, ++cnt)
    {

      // the global id
      auto part_id = cell_parts[cnt];
        
      // get the destination rank
      size_t rank = rank_owner[part_id];
      
      // sending to another block
      if ((part_id != part_ids[b]) || (world_rank != rank))
      {
        // count erasure
        erase_counter++;
        // pack the data
        tempbuf.clear();
        this_block->pack(local_id, tempbuf, false);
        // length of buffer information
        auto len = tempbuf.size();
        // destined for this rank
        if (world_rank == rank) {
          auto i = part_map.at(part_id);
          localcounts[i] += len;
        }
        // otherwise destined for another
        else
          sendcounts[part_id] += len;
      } // migrate

    } // cells

  } // blocks

  // for every populated block, we will be sending the block id and
  // the number of bytes
  for (auto & i : sendcounts)
    if (i>0)
      i+= sizeof(int_t) + sizeof(size_t);

  //------------------------------------
  // Create storage for buffers
  
  std::vector<std::vector<byte_t>> localbuf(num_parts);
  for (int_t b=0; b<num_parts; ++b)
    localbuf[b].reserve(localcounts[b]);

  std::vector<std::vector<size_t>> block_erase_ids(num_blks);
  for (int_t b=0; b<num_blks; ++b)
    block_erase_ids.reserve(block_erase_counts[b]);

  std::vector<int> senddispls(tot_parts + 1);
  senddispls[0] = 0;
  std::partial_sum(sendcounts.begin(), sendcounts.end(), &senddispls[1]);
  
  std::vector<byte_t> sendbuf(senddispls.back());

  std::fill(sendcounts.begin(), sendcounts.end(), 0);

  //------------------------------------
  // Fill buffers 

  // add the destination block ids and the number of bytes
  for (int_t b=0; b<tot_parts; ++b) {
    size_t nbytes = senddispls[b+1] - senddispls[b];
    if (nbytes) {
      auto dst = sendbuf.data() + senddispls[b];
      *reinterpret_cast<size_t*>(dst) = nbytes;
      dst += sizeof(size_t);
      *reinterpret_cast<int_t*>(dst) = b;
      sendcounts[b] += sizeof(int_t) + sizeof(size_t);
    }
  }

  // now add the rest
  for (int_t b=0, cnt=0; b<num_blks; ++b) {

    auto this_block = dynamic_cast<mesh_unstruct_block_t*>(blocks_[b]);

    size_t num_elements = this_block->local2global_.at(num_dims_).size();
    auto & erase_local_ids = block_erase_ids[b];

    for(size_t local_id(0); local_id < num_elements; ++local_id, ++cnt)
    {

      // the global id
      auto part_id = cell_parts[cnt];
      
      // get the destination rank
      size_t rank = rank_owner[part_id];
      
      // sending to another block
      if ((part_id != part_ids[b]) || (world_rank != rank))
      {
        // mark for deletion
        erase_local_ids.emplace_back(local_id);
        // get the destination rank
        size_t rank = rank_owner[part_id];
        // get the buffer the buffer
        tempbuf.clear();
        this_block->pack(local_id, tempbuf, false);
        // length of buffer information
        auto len = tempbuf.size();
        // destined for this rank
        if (world_rank == rank) {
          auto i = part_map.at(part_id);
          localbuf[i].insert(localbuf[i].end(), tempbuf.begin(), tempbuf.end());
        }
        // otherwise destined for another
        else {
          auto start = senddispls[part_id] + sendcounts[part_id];
          std::memcpy(sendbuf.data() + start, tempbuf.data(), tempbuf.size()); 
          sendcounts[part_id] += len;
        }
      } // migrate
    
    } // for entities
  } // blocks

  //----------------------------------------------------------------------------
  // Erase entities to be migrated
  //----------------------------------------------------------------------------
    
  for (int_t b=0; b<num_blks; ++b) {

    auto & erase_local_ids = block_erase_ids[b];
    if(erase_local_ids.size()) {
    
      auto this_block = dynamic_cast<mesh_unstruct_block_t*>(blocks_[b]);

      // sort them first
      std::sort(erase_local_ids.begin(), erase_local_ids.end());
      auto last = std::unique(erase_local_ids.begin(), erase_local_ids.end());
      if (last != erase_local_ids.end())
        THROW_ERROR("duplicate ids to delete");

      // erase mesh elements
      this_block->erase(erase_local_ids);

    } // erase
      
  }// blocks
  
  //----------------------------------------------------------------------------
  // Any migration setup
  //----------------------------------------------------------------------------

  for (int_t b=0; b<num_blks && num_dims_==3; ++b) {
    auto this_block = dynamic_cast<mesh_unstruct_block_t*>(blocks_[b]);
    this_block->sort_face_vertices();
  }
  
  
  //----------------------------------------------------------------------------
  // Send information.
  //----------------------------------------------------------------------------
  
  // collapse send counts, displacements
  if (sendcounts.size() < world_size) {
    sendcounts.resize(world_size);
    senddispls.resize(world_size+1);
  }

  for (size_t r=0; r<world_size; ++r) {
    senddispls[r+1] = senddispls[ part_dist[r+1] ];
    sendcounts[r] = senddispls[r+1] - senddispls[r];
  }
  sendcounts.resize(world_size);
  senddispls.resize(world_size+1);

  // the mpi data type for size_t
  std::vector<int> recvcounts(world_size, 0);
  comm_.all_to_all(sendcounts, recvcounts);
  
  // how much info will we be receiving
  std::vector<int> recvdispls(world_size + 1);
  recvdispls[0] = 0;
  std::partial_sum(recvcounts.begin(), recvcounts.end(), &recvdispls[1]);
  
  std::vector<byte_t> recvbuf(recvdispls.back());

  // some ranks might not actualy participate
  auto num_sendrecv =
    std::accumulate(sendcounts.begin(), sendcounts.end(), 0) +
    std::accumulate(recvcounts.begin(), recvcounts.end(), 0);
  sub_comm_t subcomm(comm_, num_sendrecv);

  //------------------------------------
  // Do the exchange
  
  mpi_request_t req;
  size_t comm_size = 0;

  if (num_sendrecv) {

    comm_size = subcomm.comm().size();

    // collapse all rank data
    for (size_t r=0, r2=0; r<world_size; ++r) {
      if (subcomm.is_included(r)) {
        sendcounts[r2] = sendcounts[r];
        recvcounts[r2] = recvcounts[r];
        senddispls[r2+1] = senddispls[r2] + sendcounts[r];
        recvdispls[r2+1] = recvdispls[r2] + recvcounts[r];
        r2++;
      }
    }

    sendcounts.resize(comm_size);
    recvcounts.resize(comm_size);
    senddispls.resize(comm_size+1);
    recvdispls.resize(comm_size+1);

    // now send the actual info
    req = subcomm.comm().all_to_allv(
        sendbuf,
        sendcounts,
        senddispls,
        recvbuf,
        recvcounts,
        recvdispls);

  } // has send or recive


  //----------------------------------------------------------------------------
  // create new blocks if needed, dont worry about deleting empty
  // ones yet
  //----------------------------------------------------------------------------

  std::map<int_t, int_t> block_map;
  
  // add existing
  for (int_t b=0; b<num_blks; ++b)
    block_map[ part_ids[b] ] = b;

  // add new
  for (auto b=part_dist[world_rank]; b<part_dist[world_rank+1]; ++b)
  {
    auto res = block_map.emplace(b, block_map.size());
    if (res.second) {
      auto newb = new mesh_unstruct_block_t(num_dims_);
      blocks_.emplace_back( newb );
    }
  }

  // may have changed
  num_blks = block_map.size();

  //----------------------------------------------------------------------------
  // Unpack information.
  //----------------------------------------------------------------------------
  
  //------------------------------------
  // Unpack local data
  
  for (int_t b=0; b<num_parts; ++b) {
    
    // get the block
    auto part_id = part_dist[world_rank] + b;
    auto bid = block_map.at(part_id);
    auto this_block = dynamic_cast<mesh_unstruct_block_t*>(blocks_[bid]);

    const auto & buffer = localbuf[b];
    auto nbytes = buffer.size();

    if (nbytes) {
      // get pointer
      const byte_t * ptr = &buffer[0];
      const byte_t * end = &buffer[nbytes];
      // unpack
      while (ptr<end) this_block->unpack(ptr, false);
      // make sre we are at the right spot in the buffer
      if (ptr != end) THROW_ERROR("Local unpacking mismatch");
    } // n
  }
  
  //------------------------------------
  // Unpack other ranks data

  // need exchange results now
  req.wait();

  // Add indices to primary
  for(size_t rank(0); rank < comm_size; ++rank) {

    auto start = recvdispls[rank];
    auto end = recvdispls[rank + 1];
    const auto * buffer = &recvbuf[start];
    const auto * buffer_end = &recvbuf[end];

    // Typically, one would check for equality between an iterator address
    // and the ending pointer address.  This might cause an infinate loop
    // if there is an error in unpacking.  So I think testing on byte
    // index is safer since there is no danger of iterating past the end.
    for(auto pos = start; pos < end;) {

      // capture start of buffer
      const auto * buffer_start = buffer;

      // get the number of bytes
      size_t nbytes;
      uncast(buffer, 1, &nbytes);

      // global partition id
      int_t part_id;
      uncast(buffer, 1, &part_id);
        
      // get the desination block
      auto b = block_map.at(part_id);
      auto this_block = dynamic_cast<mesh_unstruct_block_t*>(blocks_[b]);
      
      while( static_cast<size_t>(std::distance(buffer_start, buffer)) < nbytes)
        this_block->unpack(buffer, false);

      // increment pointer
      pos += std::distance(buffer_start, buffer);

    } // for

    // make sre we are at the right spot in the buffer
    if (buffer != buffer_end) THROW_ERROR("Unpacking mismatch");

  } // for
  
  //----------------------------------------------------------------------------
  // Collapse the empty blocks
  //----------------------------------------------------------------------------
  
  // move empty blocks to back
  int_t block_count = 0;

  // collapse and clean up
  for (int_t b=0; b<num_blks; ++b) {
    auto this_block = dynamic_cast<mesh_unstruct_block_t*>(blocks_[b]);
    auto it = this_block->local2global_.find(num_dims_);
    // keep it
    if (it != this_block->local2global_.end() && it->second.size()) {
      if (b!=block_count) blocks_[block_count] = blocks_[b];
      block_count++;
    }
    // delete it
    else {
      delete blocks_[b];
    }
  }

  // delete empty blocks
  blocks_.resize(block_count);
  
  //----------------------------------------------------------------------------
  // Finalize
  //----------------------------------------------------------------------------
  
  // reorder regions and sides
  for (const auto & b : blocks_) {
    auto this_block = dynamic_cast<mesh_unstruct_block_t*>(b);
    this_block->group_regions();
    this_block->group_sides();
  }
  
  // renumber blocks DONT TRUST PARMETIS
  std::vector<int_t> block_counts(world_size);
  comm_.all_gather(block_count, block_counts);
  std::partial_sum(block_counts.begin(), block_counts.end(), &block_offsets_[1]);
}

////////////////////////////////////////////////////////////////////////////////
/// create a cell connecitivity graph
////////////////////////////////////////////////////////////////////////////////
void exodus_mesh_t::graph(
    int_t min_connections,
    std::vector<size_t> & dist,
    crs_t<size_t, int_t> & grph)
{
  auto num_parts = num_blocks();
  std::vector< const crs_t<int_t, int_t> * > cells2vert(num_parts);
  std::vector< const std::vector<long_t> * > vert_local2global(num_parts);

  for (int_t b=0; b<num_parts; ++b) {
    auto this_block =
      dynamic_cast<const mesh_unstruct_block_t*>(blocks_[b]);
    cells2vert[b] = & this_block->connectivity_.at({num_dims_, 0});
    vert_local2global[b] = & this_block->local2global_.at(0);
  }
  
  prl::graph(
      comm_,
      block_offsets_,
      cells2vert,
      vert_local2global,
      min_connections,
      dist,
      grph);
}

////////////////////////////////////////////////////////////////////////////////
/// Construct the global mapping
////////////////////////////////////////////////////////////////////////////////
void exodus_mesh_t::number()
{
  number_disjoint(num_dims_);
}
    
////////////////////////////////////////////////////////////////////////////////
/// Add a halo
////////////////////////////////////////////////////////////////////////////////
void exodus_mesh_t::_build_halo(
    int_t num_ghost,
    bool with_corners,
    std::vector<comm_map_block_t> & comm_maps)
{
    
  if (num_ghost<0) return;
  
  //----------------------------------------------------------------------------
  // make the graph

  std::vector<size_t> cell_dist;
  crs_t<size_t, int_t> grph;
  int_t nconn = with_corners ? 1 : num_dims_;
  
  graph(nconn, cell_dist, grph);


  //----------------------------------------------------------------------------
  // Determine primary/ghost division
  
  auto comm_rank = comm_.rank();
  auto blkstart = block_offsets_[comm_rank];

  const auto & graph_offsets = grph.offsets;
  const auto & graph_indices = grph.indices;

  std::unordered_map<long_t, int_t> ghost_map;
  std::vector<entity_info_t> ghost_info;
  std::vector<entity_info_t> shared_info;

  auto num_parts = num_blocks();

  comm_maps.reserve(num_parts);

  for ( int_t b=0, cnt=0; b<num_parts; ++b)
  {

    // starting and ending cell indices
    auto bid = blkstart + b;
    auto cellstart = cell_dist[bid];
    auto cellend = cell_dist[bid + 1];
    auto num_elem = cellend - cellstart;
  
    auto this_block = block(b);
    auto & gids = this_block->global_ids(num_dims_);

    shared_info.clear();
    ghost_info.clear();
    ghost_map.clear();

    for(size_t lid = 0; lid < num_elem; ++lid, ++cnt)
    {

      // loop over neighbors
      for(
          auto i = graph_offsets[cnt];
          i < graph_offsets[cnt + 1];
          ++i)
      {
        // who owns the neighbor cell
        auto ghost_gid = graph_indices[i];
        int_t blk = owner(cell_dist, ghost_gid);

        // neighbor is a ghost cell
        if(blk != bid) {
          // create possible new ghost
          auto res = ghost_map.emplace(
              ghost_gid,
              num_elem + ghost_map.size());
          // if inserted
          if (res.second) gids.emplace_back(ghost_gid);
          // the new local id
          auto new_lid = res.first->second;
          // determine ghost info
          auto ghost_lid = ghost_gid - cell_dist[blk];
          auto rank = owner(block_offsets_, blk); 
          ghost_info.emplace_back(new_lid, rank, blk, ghost_lid);
          // keep track of shared ranks for the parent cell
          shared_info.emplace_back(lid, rank, blk, lid);
        }
      } // neighbors

    }
    
    // create a new comm map
    auto num_ghost = ghost_map.size();
    auto num_tot_elem = num_elem + num_ghost;
    comm_maps.emplace_back(bid, num_tot_elem, shared_info, ghost_info);

    // update number of entities
    this_block->num_all_entities(num_dims_) = num_tot_elem;


  } // blocks
  

}
} // namespace
