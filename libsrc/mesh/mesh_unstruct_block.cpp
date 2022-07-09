#include "mesh_unstruct_block.hpp"
#include "connect.hpp"

#include "mesh_unstruct_block_cpu.hpp"
#include "mesh_boundary.hpp"

#include "comm/comm_map.hpp" 

#ifdef HAVE_EXODUS
#include "io/exo_writer.hpp"
#endif
#include "io/vtk_writer.hpp"
#include "math/reorder.hpp"
#include "math/subdivide.hpp"
#include "utils/cast.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Structured mesh block
////////////////////////////////////////////////////////////////////////////////
mesh_unstruct_block_t::mesh_unstruct_block_t(int_t num_dims) : 
  mesh_block_t(num_dims)
{
  // make sure it doesnt accidentally get deleted
  desired_connectivity_.emplace(num_dims, 0);
  if (num_dims==3) {
    desired_connectivity_.emplace(num_dims-1, 0);
    desired_connectivity_.emplace(num_dims, num_dims-1);
  }

}

////////////////////////////////////////////////////////////////////////////////
/// The number of cells
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::post_connectivity()
{
  for (int_t d=1; d<num_dims_; ++d)
    set_num_entities(d);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a region
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::add_region(int_t id, int_t num)
{
  regions_.emplace_back(id);
  
  if (region_cell_offsets_.empty())
    region_cell_offsets_.emplace_back(0);
  region_cell_offsets_.emplace_back(region_cell_offsets_.back() + num);
}

////////////////////////////////////////////////////////////////////////////////
/// create a map of sorted face vertices
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::sort_face_vertices()
{
  auto it = connectivity_.find({num_dims_-1, 0});
  if (it == connectivity_.end()) return;

  auto & faces2verts = it->second;
  auto & indices = faces2verts.indices;
  auto & offsets = faces2verts.offsets;

  sorted_face_vertices_.clear();

  std::vector<int_t> sorted_vs;
  size_t num_faces = faces2verts.size();

  for (size_t f=0; f<num_faces; ++f) {
    auto start = offsets[f];
    auto end = offsets[f+1];
    sorted_vs.assign( &indices[start], &indices[end] );
    std::sort( sorted_vs.begin(), sorted_vs.end() );
    sorted_face_vertices_.emplace(sorted_vs, f);
  }
  
  assert( num_faces == sorted_face_vertices_.size() );
}
  
////////////////////////////////////////////////////////////////////////////////
/// Build extra side connectivity
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::build_cell_sides()
{
  cell_sides_.clear();
  
  auto num_sides = side_id_.size();
  const auto & local2global = local2global_.at(num_dims_);

  for (size_t s=0; s<num_sides; ++s) {
    auto lid = side_cell_[s];
    auto gid = local2global.at(lid);
    cell_sides_[gid].emplace_back(s);
  }

}

////////////////////////////////////////////////////////////////////////////////
/// extract boundaries
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::extract_boundary(
    int_t boundary_id,
    std::vector<boundary_face_t> & bfaces)
{
  // clear any existing results
  bfaces.clear();

  int_t side_id;
  auto has_boundary = find_boundary_from_id(boundary_id, side_id);
  if (!has_boundary) return;
  
  auto it = std::find(side_sets_.begin(), side_sets_.end(), side_id);
  if (it == side_sets_.end()) return;

  auto ss = std::distance(side_sets_.begin(), it);
  
  const auto & cell2faces = connectivity(num_dims_, num_dims_-1);
  const auto & c2f_offsets = cell2faces.offsets;
  const auto & c2f_indices = cell2faces.indices;
  
  auto side_start = side_set_offsets_[ss];
  auto side_end = side_set_offsets_[ss+1];
  auto num_side = side_end - side_start;
  
  bfaces.reserve(num_side);
  
  for (auto s=side_start; s<side_end; ++s) {
    auto c = side_cell_[s];
    auto f = side_cell_face_[s];
    auto offset = c2f_offsets[c];
    auto fid = c2f_indices[offset + f];
    bfaces.emplace_back(fid, 1, c);
  } // sides

}

////////////////////////////////////////////////////////////////////////////////
/// pack a cell into a buffer
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::pack(
    int_t local_id,
    std::vector<byte_t>& buffer,
    bool as_ghost) const
{
  auto & cell_local2global_ = local2global_.at(num_dims_);
  auto & vertex_local2global_ = local2global_.at(0);
  
  //--------------------------------------------------------------------------
  // pack cell info
  
  // add mesh global id to buffer
  long_t global_id = cell_local2global_[local_id];
  cast_insert( &global_id, 1, buffer );
  // add cell block info
  cast_insert( &cell_region_[local_id], 1, buffer );
  cast_insert( &cell_type_[local_id], 1, buffer );
  
  //--------------------------------------------------------------------------
  // pack vertices
  {
    auto & cell_vertices_ = connectivity_.at({num_dims_, 0});
    auto & offsets = cell_vertices_.offsets;
    auto & indices = cell_vertices_.indices;
  
    // get number of vertices
    const auto start = offsets[local_id];
    const auto end   = offsets[local_id+1];
    int_t num_verts = end - start;
    // add num verts to buffer
    cast_insert( &num_verts, 1, buffer );
    // add global vertex indices to buffer
    for (auto i=start; i<end; ++i) {
      auto lid = indices[i];
      long_t gid = vertex_local2global_[lid];
      cast_insert( &gid, 1, buffer );
    }
    // also pack coordinates too, not ideal but easiest
    for (auto i=start; i<end; ++i) {
      auto lid = indices[i];
      const real_t * coord = vertex_coords_.data() + lid*num_dims_;
      cast_insert( coord, num_dims_, buffer );
    }
  }
  
  //--------------------------------------------------------------------------
  // pack faces
  if (num_dims_==3) {
    
    auto & cell_faces_ = connectivity_.at({3, 2});
    auto & c2f_offsets = cell_faces_.offsets;
    auto & c2f_indices = cell_faces_.indices;
    
    auto & face_vertices_ = connectivity_.at({2, 0});
    auto & f2v_offsets = face_vertices_.offsets;
    auto & f2v_indices = face_vertices_.indices;

    std::vector<size_t> vs;
      
    // get number of faces
    const auto f_start = c2f_offsets[local_id];
    const auto f_end = c2f_offsets[local_id + 1];
    int_t num_faces = f_end - f_start;
    // add num faces to buffer
    cast_insert( &num_faces, 1, buffer );
    // loop over faces
    for ( auto f=f_start; f<f_end; ++f ) {
      auto fid = c2f_indices[f];
      // get face verts
      const auto v_start = f2v_offsets[fid];
      const auto v_end = f2v_offsets[fid+1];
      auto num_verts = v_end - v_start;
      vs.assign( &f2v_indices[v_start], &f2v_indices[v_end] );
      // if this cell is not the owner, then reverse the face verts
      if ( face_owner_[fid] != local_id )
        std::reverse( vs.begin(), vs.end() );
      // add number of verts to buffer
      cast_insert( &num_verts, 1, buffer );
      // add face verts to buffer
      for ( auto vid : vs ) {
        auto gid = vertex_local2global_[vid];
        cast_insert( &gid, 1, buffer );
      }
    } // face

  } // dims == 3
  
  //--------------------------------------------------------------------------
  // pack side sets

  // add side_sets to buffer
  auto sit = cell_sides_.find(global_id);
  if ( sit != cell_sides_.end() && !as_ghost ) {
    const auto & sides = sit->second;
    int_t num_sides = sides.size();
    cast_insert( &num_sides, 1, buffer );
    for ( auto s : sides ) {
      cast_insert( &side_id_[s], 1, buffer );
      cast_insert( &side_cell_face_[s], 1, buffer );
    }
  }
  else {
    int_t num_sides = 0;
    cast_insert( &num_sides, 1, buffer );
  } // has sides
}

////////////////////////////////////////////////////////////////////////////////
/// erase a list of cells
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::erase(const std::vector<size_t> & local_ids)
{
  if (local_ids.empty()) return;

  // assume sorted
  assert( std::is_sorted( local_ids.begin(), local_ids.end() )
      && "entries to delete are not sorted" );

  //--------------------------------------------------------------------------
  // Erase elements
  
  auto & cell_local2global_ = local2global_.at(num_dims_);
  auto & cell_global2local_ = global2local_.at(num_dims_);

  // erase any sides that are no longer used
  std::vector<int_t> delete_sides;
  
  auto num_cells = cell_local2global_.size();
  size_t num_remove = 0;
  
  // keep track of cell renumbering
  std::vector<int_t> cell_old2new(num_cells, -1);

  auto delete_it = local_ids.begin();
  
  for (
      size_t old_local_id=0, new_local_id=0, region_id=0;
      region_id<regions_.size();
      ++region_id )
  {

    size_t region_end = region_cell_offsets_[region_id+1];
    
    for (; old_local_id < region_end; ++old_local_id)
    {
      
      auto global_id = cell_local2global_[old_local_id];
      
      // skip deleted items
      if ( delete_it != local_ids.end() ) {
        if ( *delete_it == old_local_id ) {
      
          cell_global2local_.erase( global_id );

          auto it = cell_sides_.find( global_id );
          if ( it != cell_sides_.end() ) {
            delete_sides.insert(
                delete_sides.end(),
                it->second.begin(),
                it->second.end() );
            cell_sides_.erase( it );
          }

          delete_it++;
          num_remove++;
          continue;
        }
      }

      // keep otherwise
      cell_local2global_[new_local_id] = cell_local2global_[old_local_id];
      cell_global2local_[global_id] = new_local_id;
      //cell_global_ids_[new_local_id] = cell_global_ids_[old_local_id];
      cell_type_[new_local_id] = cell_type_[old_local_id];
      cell_region_[new_local_id] = cell_region_[old_local_id];
      cell_old2new[ old_local_id ] = new_local_id;
      new_local_id ++;

    } // cells

    // update cell region offsets
    region_cell_offsets_[region_id+1] = new_local_id;

  } // regions

  // resize
  num_cells -= num_remove;
  cell_local2global_.resize(num_cells);
  cell_type_.resize(num_cells);
  cell_region_.resize(num_cells);
  
  // collapse regions
  size_t num_regions=0;
  for (size_t r=0; r<regions_.size(); ++r) {
    auto n = region_cell_offsets_[r+1] - region_cell_offsets_[r]; 
    if (n) {
      region_cell_offsets_[num_regions+1] = region_cell_offsets_[num_regions] + n;
      regions_[num_regions] = regions_[r];
      num_regions++;
    }
  }
  region_cell_offsets_.resize(num_regions+1);
  regions_.resize(num_regions);

  // erase cells to vertices info
  auto & cell_vertices_ = connectivity_.at({num_dims_, 0});
  cell_vertices_.erase(local_ids);

  //--------------------------------------------------------------------------
  // Erase faces
  
  // temporary storage 
  std::vector<int_t> counts;
  std::vector<int_t> old2new;
    
  crs_t<int_t, int_t> * face_vertices_ = nullptr;

  if (num_dims_ == 3) {
    
    // erace cell to face info
    auto & cell_faces_ = connectivity_.at({3, 2});
    cell_faces_.erase(local_ids);
    
    auto & c2f_offsets = cell_faces_.offsets;
    auto & c2f_indices = cell_faces_.indices;
    
    face_vertices_ = &connectivity_.at({2, 0});
    auto & f2v_offsets = face_vertices_->offsets;
    auto & f2v_indices = face_vertices_->indices;

    //-----------------------------------
    // Determine unused faces
    size_t num_faces = face_vertices_->size();
    counts.resize( num_faces, 0 );

    for ( auto i : c2f_indices ) counts[i]++; 

    bool has_unused_faces = false;
    for ( decltype(num_faces) i=0; i<num_faces; ++i ) {
      if ( counts[i] == 0 ) {
        has_unused_faces = true;
        break;
      }
    }
    
    //-----------------------------------
    // Delete faces
    
    if ( has_unused_faces ) {
      
      // storage for renumberings
      old2new.resize( num_faces, -1 );

      std::vector<int_t> deleted_faces;
      deleted_faces.reserve( num_faces );

      for ( size_t old_local_id=0, new_local_id=0; old_local_id<num_faces; ++old_local_id )
      {
        
        // skip deleted items
        if ( counts[old_local_id] == 0 ) {
            deleted_faces.emplace_back( old_local_id );
            continue;
        }

        // keep otherwise
        face_owner_[new_local_id] = face_owner_[old_local_id];

        old2new[ old_local_id ] = new_local_id;
        
        new_local_id++;
        
      }
   
      // resize
      num_faces -= deleted_faces.size();
      face_owner_.resize(num_faces);
   
      // delete face to vertex connectivity
      face_vertices_->erase(deleted_faces);

      // renumber cell connectivity
      for ( auto & i : c2f_indices ) {
        i = old2new[i];
        assert( i<static_cast<int_t>(num_faces) && "Messed up renumbering #1" );
      }


      // reverse any vertices whose owner was deleted
      for ( size_t f=0; f<num_faces; ++f ) {
        // get the old owner
        auto owner_id = face_owner_[f];
        auto new_owner_id = cell_old2new[owner_id];
        // if it isn't the global list, reverse it
        if ( new_owner_id == -1 ) {
          auto start = f2v_offsets[f];
          auto end = f2v_offsets[f+1];
          std::reverse( &f2v_indices[start], &f2v_indices[end] );
          // reset the id
          face_owner_[f] = -1;
        }
        // otherwise, renumber it
        else {
          face_owner_[f] = new_owner_id;
        }
      }
   
      // figure out any of the reset face owners
      for ( size_t c=0; c<num_cells; ++c ){
        auto start = c2f_offsets[c];
        auto end = c2f_offsets[c+1];
        for ( auto f=start; f<end; ++f ) {
          auto fid = c2f_indices[f];
          if ( face_owner_[fid] == -1 )
            face_owner_[fid] = c;
        }
      }
      

    } // has unused faces

  } // dims == 3
  
  //--------------------------------------------------------------------------
  // Erase vertices

  //-------------------------------------
  // Determine unused vertices
  
  int_t num_vertices = vertex_coords_.size() / num_dims_;
  counts.clear();
  counts.resize( num_vertices, 0 );

  for ( auto i : cell_vertices_.indices ) counts[i]++; 

  bool has_unused_vertices = false;
  for ( int_t i=0; i<num_vertices; ++i ) {
    if ( counts[i] == 0 ) {
      has_unused_vertices = true;
      break;
    }
  }
  
  //-------------------------------------
  // Delete vertices
  
  if ( has_unused_vertices ) {
  
    auto & vertex_local2global_ = local2global_.at(0);
    auto & vertex_global2local_ = global2local_.at(0);
  
    // storage for renumberings
    old2new.clear();
    old2new.resize( num_vertices, -1 );
  
    size_t deleted_vertices = 0;

    for (
        int_t old_local_id=0, new_local_id=0;
        old_local_id<num_vertices;
        ++old_local_id )
    {
  
      // get global and local id
      auto global_id = vertex_local2global_[old_local_id];
    
      // skip deleted items
      if ( counts[old_local_id] == 0 ) {
          vertex_global2local_.erase( global_id );
          deleted_vertices++;
          continue;
      }

      // keep otherwise
      for ( int d=0; d<num_dims_; ++d )
        vertex_coords_[ new_local_id*num_dims_ + d ] =
          vertex_coords_[ old_local_id*num_dims_ + d ];
      vertex_local2global_[new_local_id] = vertex_local2global_[old_local_id];
      vertex_global2local_[global_id] = new_local_id;

      old2new[ old_local_id ] = new_local_id;
      
      new_local_id++;
      
    }

    // resize
    int_t vertex_count = num_vertices - deleted_vertices;
    vertex_local2global_.resize(vertex_count);
    vertex_coords_.resize(vertex_count*num_dims_);
    
    // renumber cell connectivity
    for ( auto & i : cell_vertices_.indices ) {
      i = old2new[i];
      assert( i > -1 && i < vertex_count && "Messed up renumbering #1" );
    }
    // renumber face connectivity
    if (face_vertices_) {
      for ( auto & i : face_vertices_->indices ) {
        i = old2new[i];
        assert( i > -1 && i < num_vertices && "Messed up renumbering #2" );
      }
    }

  } // delete vertices

    
  // delete the sorted face vertices.  you cant change the key of a map, so
  // you cant renumber the vertices.  just blow away the whole map and
  // rebuild it.
  if (num_dims_ == 3) sort_face_vertices();
  
  //--------------------------------------------------------------------------
  // Delete sides

  if ( !delete_sides.empty() ) {

    std::sort( delete_sides.begin(), delete_sides.end() );

    size_t deleted_sides = 0;

    // storage for renumberings
    int_t num_sides = side_id_.size();
    std::vector<int_t> old2new( num_sides, -1 );
    
    auto delete_it = delete_sides.begin();
    for (
        int_t old_local_id=0, new_local_id=0;
        old_local_id<num_sides;
        ++old_local_id )
    {

      // skip deleted items
      if ( delete_it != delete_sides.end() ) {
        if ( *delete_it == old_local_id ) {
          ++deleted_sides;
          ++delete_it;
          continue;
        }
      }

      // keep otherwise
      side_id_[new_local_id] = side_id_[old_local_id];
      side_cell_[new_local_id] = side_cell_[old_local_id];
      side_cell_face_[new_local_id] = side_cell_face_[old_local_id];
      old2new[ old_local_id ] = new_local_id;
      new_local_id++;

    }

    // can resize side ids now
    num_sides -= deleted_sides;
    side_id_.resize( num_sides );
    side_cell_.resize( num_sides );
    side_cell_face_.resize( num_sides );

    // renumber element sides
    for (auto & e : side_cell_)
      e = cell_old2new[e];

    // renumber element sides
    for ( auto & pair : cell_sides_ )
      for ( auto & s : pair.second )
        s = old2new[s];

  }

}

////////////////////////////////////////////////////////////////////////////////
/// unpack a cell
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::unpack(const byte_t * & buffer, bool as_ghost)
{
  auto & cell_local2global_ = local2global_[num_dims_];
  auto & cell_global2local_ = global2local_[num_dims_];
  
  auto & vertex_local2global_ = local2global_[0];
  auto & vertex_global2local_ = global2local_[0];
  
  // compute the new local id
  int_t local_id = cell_local2global_.size();
  // get global mesh id
  long_t global_id;
  uncast( buffer, 1, &global_id );
  // add mapping
  cell_local2global_.emplace_back( global_id );
  auto res = cell_global2local_.emplace( std::make_pair(global_id, local_id) );
  assert( res.second && "Global id already exists but shouldn't" );
  UNUSED(res); // make production build happy
  
  // get cell info
  int_t cell_region;
  uncast( buffer, 1, &cell_region );
  cell_region_.emplace_back(cell_region);

  shape_t cell_type;
  uncast( buffer, 1, &cell_type );
  cell_type_.emplace_back(cell_type);
  
  //--------------------------------------------------------------------------
  // Vertices

  // get the number of vertices
  int_t num_verts;
  uncast( buffer, 1, &num_verts );
  
  // retrieve global vertex ids
  std::vector<long_t> vs(num_verts);
  uncast( buffer, num_verts, vs.data() );
  
  // unpack vertex coordinates
  std::vector<real_t> coords(num_dims_*num_verts);
  uncast( buffer, coords.size(), coords.data() );
  
  // need to convert global vertex ids to local ids
  vertex_local2global_.reserve( vertex_local2global_.size() + num_verts);
  vertex_coords_.reserve( vertex_coords_.size() + num_verts*num_dims_ );

  std::vector<int_t> new_vs;
  new_vs.reserve(num_verts);

  for ( int_t i=0; i<num_verts; ++i ) {
    auto gid = vs[i];
    auto lid = vertex_global2local_.size();
    auto it = vertex_global2local_.find(gid);
    // doesnt exist, so add entries to maps
    if (it == vertex_global2local_.end()) {
      if (!as_ghost) {
        vertex_global2local_.emplace(gid, lid);
        vertex_local2global_.emplace_back( gid );
        for ( int_t dim=0; dim<num_dims_; ++dim )
          vertex_coords_.emplace_back( coords[i*num_dims_ + dim] );
        new_vs.emplace_back(lid);
      }
    }
    // exists, so get local id
    else {
      new_vs.emplace_back(it->second);
    }
  }

  // add vertex connectivity
  auto & cell_vertices_ = connectivity_[{num_dims_, 0}];
  cell_vertices_.push_back(new_vs.begin(), new_vs.end());

  //--------------------------------------------------------------------------
  // faces
  if (num_dims_ == 3) {
    
    auto & cell_faces_ = connectivity_[{3, 2}];
    auto & face_vertices_ = connectivity_[{2, 0}];

    // get the number of faces
    int_t num_faces;
    uncast( buffer, 1, &num_faces );

    // storage to keep track of new cell faces
    std::vector<size_t> fs;
    fs.reserve( num_faces );
      
    // storage for vertices
    std::vector<size_t> vs;
    std::vector<int_t> sorted_vs;

    // loop over faces and unpack verts
    for ( int_t i=0; i<num_faces; ++i ) {

      // the number of vertices
      int_t num_verts;
      uncast( buffer, 1, &num_verts );

      // the ids of the vertices
      vs.clear();
      vs.resize(num_verts);
      uncast( buffer, num_verts, vs.data() );
  
      // convert vertex global ids to local ids
      bool has_ghost = false;
      for ( auto & v : vs ) {
        auto it = vertex_global2local_.find(v);
        if (it != vertex_global2local_.end())
          v = it->second;
        else {
          has_ghost = true;
          break;
        }
      }

      // skip face if outside domain
      if (has_ghost) continue;
      
      // sort vertices for search
      sorted_vs.assign(vs.begin(), vs.end());
      std::sort( sorted_vs.begin(), sorted_vs.end() );
      
      // search for the face to see if it exists
      auto it = sorted_face_vertices_.find( sorted_vs );
      // no match, add it to the connectivity list
      if ( it == sorted_face_vertices_.end() ) {
        if (!as_ghost) {
          auto new_face_id = face_vertices_.size();
          fs.emplace_back( new_face_id );
          face_vertices_.push_back( vs.begin(), vs.end() );
          face_owner_.emplace_back( local_id );
          sorted_face_vertices_.emplace( sorted_vs, new_face_id );
          assert( static_cast<size_t>(face_vertices_.size()) == sorted_face_vertices_.size() );
        }
      }
      // if there was a match, use the id
      else {
        fs.emplace_back(it->second);
      }

    }
    
    // now append to the connectivity list
    cell_faces_.push_back( fs.begin(), fs.end() );

  } // faces

  //--------------------------------------------------------------------------
  // Side sets

  int_t num_sides;
  uncast( buffer, 1, &num_sides );
    
  if ( num_sides > 0 ) {
    auto & sides = cell_sides_[global_id];
    sides.reserve(sides.size() + num_sides);
    for ( int_t i=0; i<num_sides; ++i ) {
      side_cell_.emplace_back(local_id);
      int_t side_id;
      uncast( buffer, 1, &side_id );
      sides.emplace_back( side_id_.size() );
      side_id_.emplace_back( side_id );
      int_t side_face;
      uncast( buffer, 1, &side_face );
      side_cell_face_.emplace_back( side_face );
    }
  } // has sides

}

////////////////////////////////////////////////////////////////////////////////
/// reorder cells to grooup regions together
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::group_regions()
{
  
  // count regions
  std::map<int_t, int_t> region_counts;
  for (auto r : cell_region_)
    region_counts[r]++;
  
  // map existing regions
  std::map<int_t, int_t> region_map;
  for (size_t r=0, rnew=0; r<regions_.size(); ++r) {
    auto rid = regions_[r];
    if (region_counts.count(rid)) {
      region_map.emplace(rid, rnew);
      rnew++;
    }
  }

  // add new regions, new ones should be in sorted order
  for (auto & rc : region_counts)
    region_map.emplace(rc.first, region_map.size());
  
  // now sort by region id
  auto num_cells = cell_region_.size();
  std::vector<int_t> order(num_cells);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(
      order.begin(),
      order.end(), 
      [&](auto a, auto b)
      { return region_map.at(cell_region_[a]) < region_map.at(cell_region_[b]); }
  );

  // now reorder
  if (!std::is_sorted(order.begin(), order.end()))
    reorder_cells(order);

  // set the new region meta data
  int_t num_regions = region_counts.size();
  regions_.resize(num_regions);
  region_cell_offsets_.resize(num_regions+1);
  
  for (const auto & rc : region_counts) {
    auto rid = rc.first;
    auto count = rc.second;
    auto r = region_map.at(rid);
    regions_[r] = rid;
    region_cell_offsets_[r+1] = count;
  }

  std::partial_sum(
      region_cell_offsets_.begin(),
      region_cell_offsets_.end(),
      region_cell_offsets_.begin());

}

////////////////////////////////////////////////////////////////////////////////
/// reorder sides together
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::group_sides()
{
  
  // count regions
  std::map<int_t, int_t> side_counts;
  for (auto s : side_id_)
    side_counts[s]++;

  int_t num_sets = side_counts.size();
  
  // now sort by side id
  auto num_sides = side_id_.size();
  std::vector<int_t> order(num_sides);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(
      order.begin(),
      order.end(), 
      [&](auto a, auto b)
      { return side_id_[a] < side_id_[b]; }
  );

  // now reorder
  if (!std::is_sorted(order.begin(), order.end()))
    reorder_sides(order);

  // set the new side meta data
  side_sets_.resize(num_sets);
  side_set_offsets_.resize(num_sets+1);
  side_set_offsets_[0] = 0;
  
  int_t s = 0;
  for (auto & pair : side_counts) {
    auto count = pair.second;
    auto ss_id = pair.first;
    side_sets_[s] = ss_id;
    side_set_offsets_[s+1] = side_set_offsets_[s] + count;
    ++s;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// reorder cells to grooup regions together
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::reorder_cells(const std::vector<int_t> & order)
{

  // determine the mapping
  std::vector<int_t> mapping(order.size());
  for (size_t i=0; i<order.size(); ++i)
    mapping[order[i]] = i;

  // reorder arrays
  auto order_start = order.begin();
  
  reorder(cell_type_, order_start);
  reorder(cell_region_, order_start);

  auto & local2global = local2global_.at(num_dims_);
  reorder(local2global, order_start);
  
  auto & global2local = global2local_.at(num_dims_);
  for (auto & pair : global2local)
    pair.second = mapping[pair.second];

  // reorder all connectivity from cells to entities
  for (int_t to_dim=0; to_dim<num_dims_; ++to_dim) {
    auto it = connectivity_.find({num_dims_, to_dim});
    if (it != connectivity_.end())
      it->second.reorder(order_start);
  }
  
  // reorder all connectivity from entities to cells
  for (int_t from_dim=0; from_dim<num_dims_; ++from_dim) {
    auto it = connectivity_.find({from_dim, num_dims_});
    if (it != connectivity_.end()) {
      for (auto & i : it->second.indices)
        i = mapping[i];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
/// reorder sides to grooup sets together
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::reorder_sides(const std::vector<int_t> & order)
{

  // determine the mapping
  std::vector<int_t> mapping(order.size());
  for (size_t i=0; i<order.size(); ++i)
    mapping[order[i]] = i;

  // reorder arrays
  auto order_start = order.begin();
  
  reorder(side_id_, order_start);
  reorder(side_cell_, order_start);
  reorder(side_cell_face_, order_start);
  
  for (auto & pair : cell_sides_)      
    for (auto & s : pair.second)
      s = mapping[s];

}


////////////////////////////////////////////////////////////////////////////////
/// compute cell midpoints
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::cell_midpoints(std::vector<real_t> & cx)
{
  auto & cell_vertices_ = connectivity_.at({num_dims_, 0});
  auto num_cells = cell_vertices_.size();
  cx.reserve(num_dims_ * num_cells);

  const auto & offsets = cell_vertices_.offsets;
  const auto & indices = cell_vertices_.indices;

  std::vector<real_t> tmp(num_dims_);

  for (int_t c=0; c<num_cells; ++c) {

    // zero temporary
    for (int_t d=0; d<num_dims_; ++d) tmp[d] = 0;

    // sum vertex coords
    auto start = offsets[c];
    auto end = offsets[c+1];
    for (auto v=start; v<end; ++v) {
      auto vid = indices[v];
      for (int_t d=0; d<num_dims_; ++d)
        tmp[d] += vertex_coords_[vid*num_dims_ + d];
    }
    
    // average result
    auto nverts = end - start;
    auto fact = static_cast<real_t>(1) / nverts;
    for (int_t d=0; d<num_dims_; ++d)
      cx.emplace_back( tmp[d] * fact );

  }
}

////////////////////////////////////////////////////////////////////////////////
/// build connecitivy
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::build_neighbors()
{

  // build whatever neighbor info is requested
  for (int_t dim = num_dims_+1; dim-- > 0;) {
    for (int_t thru=0; thru<=num_dims_; ++thru) {
      if (thru==dim) continue;
      build_neighbors(dim, thru);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// compute neighbors
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::build_neighbors(int_t dim, int_t thru) 
{
  auto mx = std::max(dim, thru);
  auto mn = std::min(dim, thru);

  build_connectivity(mx, mn);
  build_connectivity(mn, mx);

  const auto & conn = connectivity(dim, thru);
  const auto & rconn = connectivity(thru, dim);

  auto & neigh = neighbors_[{dim, thru}];
  prl::neighbors(conn, rconn, neigh);
}



////////////////////////////////////////////////////////////////////////////////
/// build connecitivy
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::build_connectivity()
{
  
  // build whatever is requested
  for (int_t from = num_dims_+1; from-- > 0;) {
    for (int_t to=0; to<=num_dims_; ++to) {
      if (to==from) continue;

      auto entry = std::make_pair(from, to);
      if (desired_connectivity_.count(entry))
        if (connectivity_.count(entry) == 0)
          build_connectivity(from, to);
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
/// compute connectivity
/// \remark This is a modified version of Logg's algorithm to account
///         for polyhedra
/// \see Logg, "Efficient representation of computational meshes",
///      Int J. Computational Science and Engineering, Vol 4, No 4, 2009.
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::build_connectivity(int_t from, int_t to) 
{
  // Make sure entities->vertices constructed constructed
  if (from>0) {
    auto & from_to_verts = connectivity_[{from, 0}];
    if (from_to_verts.empty()) build(from);
  }

  if (to>0) {
    auto & to_to_verts = connectivity_[{to, 0}];
    if (to_to_verts.empty()) build(to);
  }

  // Return if connectivity is already constructed
  auto & conn = connectivity_[{from, to}];
  if (conn.size()) return;

  // if from is less than to, flip
  if (from < to) {
    build_connectivity(to, from);
    transpose(from, to);
  }

  // otherwise compute connectivity
  else {
    int_t through = (from==0 && to==0) ? num_dims_ : 0;
    build_connectivity(from, through);
    build_connectivity(through, to);
    intersect(from, to, through);
  }

}

////////////////////////////////////////////////////////////////////////////////
/// Build cell -> entity, entity->vertex connectivity for 0<d<num_dims
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::transpose(int_t from, int_t to) 
{
  if (from>=to) return;
  prl::transpose( connectivity_[{to, from}], connectivity_[{from, to}]);
}

////////////////////////////////////////////////////////////////////////////////
/// Build cell -> entity, entity->vertex connectivity for 0<d<num_dims
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::build(int_t dim) 
{
  if (dim<=0 || dim>=num_dims_) return;

  auto & ents2vertices = connectivity_[{dim,0}];
  auto & cells2ents = connectivity_[{num_dims_, dim}];
  auto & cell2vertices = connectivity_.at({num_dims_, 0});
  std::map<std::vector<int_t>, int_t> sorted_vertices_to_ents;

  int_t cell_counter = 0;
  auto num_owned = num_owned_cells();

  // build the connecitivity array
  prl::build(
      cell2vertices, cells2ents, ents2vertices,
      sorted_vertices_to_ents,
      
      [=,&cell_counter](const auto & vs, auto & ent_vs) {

        //--- If three or more, this must be a closed poly (hopefully)
        if (vs.size() >= 3) {
          for (
            auto v0 = vs.begin(), v1 = std::next(v0);
            v0 != vs.end();
            ++v0, ++v1
          ) {
            if (v1 == vs.end()) { v1 = vs.begin(); }
            ent_vs.push_back({*v0, *v1});
          }
        }
        //---  Ghosts might have an insufficient number of vertices
        else if (vs.size() > static_cast<size_t>(dim)) {
          ent_vs.push_back(vs);
        }
        //--- Done
        
        cell_counter++;

        return cell_counter <= num_owned;

      });
}

////////////////////////////////////////////////////////////////////////////////
/// Compute d->d' from d->d'' and d''->d' for d>=d
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::intersect(int_t from, int_t to, int_t through) 
{
  if (from<to) return;
  auto & from2to = connectivity_[{from, to}];
  const auto & from2through = connectivity_[{from, through}];
  const auto & through2to = connectivity_[{through, to}];
  prl::intersect(from2through, through2to, from2to);
}

////////////////////////////////////////////////////////////////////////////////
/// request geeometry
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::request_face_geometry()
{
  if (num_dims_ == 1) {
    request_connectivity(1, 0);
    request_connectivity(0, 1);
  }
  else if (num_dims_ == 2) {
    request_connectivity(1, 0);
  }
  else {
    request_connectivity(2, 0);
  }
}

void mesh_unstruct_block_t::request_cell_geometry()
{
  if (num_dims_ == 1) {
    request_connectivity(1, 0);
  }
  else if (num_dims_ == 2) {
    request_connectivity(2, 0);
  }
  else {
    request_connectivity(2, 0);
    request_connectivity(3, 2);
  }
}

////////////////////////////////////////////////////////////////////////////////
/// builld geeometry
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::_build_geometry(bool with_face_geom, bool with_cell_geom)
{

  auto num_owned = num_owned_cells();

  if (num_dims_ == 1) {
    const auto & cell2verts = connectivity(1, 0);
    const auto & vert2cells = connectivity(0, 1);
    if (with_face_geom)
      unstruct_face_geometry_1d(
          vert2cells.size(),
          num_owned,
          cell2verts.offsets.data(),
          cell2verts.indices.data(),
          vert2cells.offsets.data(),
          vert2cells.indices.data(),
          vertex_coords_.data(),
          face_centroids_.data(),
          face_normals_.data(),
          face_areas_.data());
    if (with_cell_geom)
       unstruct_cell_geometry_1d(
           num_owned,
           cell2verts.offsets.data(),
           cell2verts.indices.data(),
           vertex_coords_.data(),
           cell_centroids_.data(),
           cell_volumes_.data());
  }
  else if (num_dims_ == 2) {
    const auto & face2verts = connectivity(1, 0);
    const auto & cell2verts = connectivity(2, 0);
    if (with_face_geom)
      unstruct_face_geometry_2d(
          face2verts.size(),
          face2verts.offsets.data(),
          face2verts.indices.data(),
          vertex_coords_.data(),
          face_centroids_.data(),
          face_normals_.data(),
          face_areas_.data());
    if (with_cell_geom)
      unstruct_cell_geometry_2d(
          num_owned,
          cell2verts.offsets.data(),
          cell2verts.indices.data(),
          vertex_coords_.data(),
          cell_centroids_.data(),
          cell_volumes_.data());
  }
  else {
    const auto & face2verts = connectivity(2, 0);
    const auto & cell2faces = connectivity(3, 2);
    if (with_face_geom)
      unstruct_face_geometry_3d(
          face2verts.size(),
          face2verts.offsets.data(),
          face2verts.indices.data(),
          vertex_coords_.data(),
          face_centroids_.data(),
          face_normals_.data(),
          face_areas_.data());
    if (with_cell_geom)
      unstruct_cell_geometry_3d(
          num_owned,
          cell2faces.offsets.data(),
          cell2faces.indices.data(),
          face2verts.offsets.data(),
          face2verts.indices.data(),
          face_owner_.data(),
          vertex_coords_.data(),
          cell_centroids_.data(),
          cell_volumes_.data());
  } // ndims
}

////////////////////////////////////////////////////////////////////////////////
/// vtk output
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::output(vtk_writer_t & vtk)  const
{
  vtk.write_points(vertex_coords_, num_dims_);
  
  auto ncells = num_owned_cells();

  vtk_cell_type_t ty;
  if (num_dims_ == 1)
    ty = vtk_cell_type_t::line;
  else if (num_dims_ == 2)
    ty = vtk_cell_type_t::polygon;
  else
    ty = vtk_cell_type_t::polyhedron;
  
  std::vector<short> cell_types(ncells, static_cast<std::uint8_t>(ty));
    
  auto & cell2verts = connectivity_.at({num_dims_, 0});
  auto c2v_offsets = make_array_ref(cell2verts.offsets.data(), ncells+1);
  auto c2v_indices = make_array_ref(cell2verts.indices.data(), c2v_offsets.back());
  
  if (num_dims_ < 3) {
    vtk.write_elements(c2v_offsets, c2v_indices, cell_types);
  }
  else {
    auto & cell2faces = connectivity_.at({num_dims_, num_dims_-1});
    auto & face2verts = connectivity_.at({num_dims_-1, 0});
  
    auto c2f_offsets = make_array_ref(cell2faces.offsets.data(), ncells+1);
    auto c2f_indices = make_array_ref(cell2faces.indices.data(), c2f_offsets.back());

    vtk.write_elements(
        c2f_offsets,
        c2f_indices,
        c2v_offsets,
        c2v_indices,
        face2verts.offsets,
        face2verts.indices,
        cell_types);
  }


}

////////////////////////////////////////////////////////////////////////////////
/// exo output
////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_EXODUS
void mesh_unstruct_block_t::output(
    exo_writer_t & exo,
    const std::string & label,
    const std::vector<side_set_info_t> & side_sets)  const
{

  auto num_cells = num_owned_cells();
  if (num_cells==0) return;

  auto num_verts = num_owned_vertices();
  int_t num_regions = regions_.size();
  auto num_sets = side_sets.size();

  const auto & cell2vert = connectivity_.at({num_dims_, 0});
  const auto & cell_local2global = local2global_.at(num_dims_);
  const auto & vert_local2global = local2global_.at(0);
  
  // write param:
  auto params = exo_writer_t::make_params(label.c_str());
  params.num_dim = num_dims_;
  params.num_nodes = num_verts;
  params.num_elem_blk = num_regions;
  params.num_elem = num_cells; 
  params.num_side_sets = num_sets;
  exo.write_params(params);
      
  // check the integer type used in the exodus file
  auto int64 = exo.is_int64();
  
  //----------------------------------------------------------------------------
  // Write each cell region
      
  for (int_t r=0; r<num_regions; ++r) {
      
    // renumber cells + verts    
    auto region_start = region_cell_offsets_[r];
    auto num_region_cells = region_cell_offsets_[r+1] - region_start;

    // helper function for getting cell vertices
    auto cell_vertices_func = [&](auto c, auto & vert_list) {
      auto local_id = region_start + c;
      const auto & vs = cell2vert.at(local_id);
      vert_list.insert(vert_list.end(), vs.begin(), vs.end());
    };
  
    // add the whole element block
    std::stringstream ss;
    ss << "Block " << id();
  
    auto type_str = exo_writer_t::geom_type(cell_type_[region_start]);
    auto region_id = regions_[r] + 1;

    if (int64) {
      exo.write_element_block<long long>(
          region_id,
          ss.str().c_str(),
          type_str,
          num_region_cells,
          cell_vertices_func);
    }
    else {
      exo.write_element_block<int>(
          region_id,
          ss.str().c_str(),
          type_str,
          num_region_cells,
          cell_vertices_func);
    }

  } // regions

  //----------------------------------------------------------------------------
  // write final element mapping
  if (int64)
    exo.write_element_map<long long>(cell_local2global);
  else
    exo.write_element_map<int>(cell_local2global);
      
  //----------------------------------------------------------------------------
  // write coordinates
  std::vector<real_t> coordinates(num_dims_ * num_verts);

  for (int_t v=0; v<num_verts; ++v) {
    for (int_t d=0; d<num_dims_; ++d )
      coordinates[v + d*num_verts] = vertex_coords_[v*num_dims_ + d];
  }
  exo.write_point_coords(coordinates);
  
  if (int64)
    exo.write_node_map<long long>(vert_local2global);
  else
    exo.write_node_map<int>(vert_local2global);
      
  //--------------------------------------------------------------------------
  // Write side sets
    
  if (num_sets) {

    std::vector<std::string> names;
    names.reserve(num_sets);

    for (const auto & side : side_sets) {

      names.emplace_back(side.label);

      auto ss_id = side.id + 1;
      auto it = std::find(side_sets_.begin(), side_sets_.end(), side.id);
      
      //--- write empty side
      if (it == side_sets_.end()) {
        exo.write_side_set_param(ss_id, 0);
      }
      //--- write side set
      else {

        auto s = std::distance(side_sets_.begin(), it);
        auto side_start = side_set_offsets_[s];
        auto num_sides = side_set_offsets_[s+1] - side_start;
    
        auto side_func = [&](auto s, auto & elem, auto & side) {
          auto local_id = side_start + s;
          elem = side_cell_[local_id];
          side = side_cell_face_[local_id];
        };
   
        if (int64)
          exo.write_side_set<long long>(ss_id, num_sides, side_func);
        else
          exo.write_side_set<int>(ss_id, num_sides, side_func);
      }

    } // id

    exo.write_side_set_names(names);

  } // sides
      
}
#endif

////////////////////////////////////////////////////////////////////////////////
/// find a cell
////////////////////////////////////////////////////////////////////////////////
bool mesh_unstruct_block_t::find_cell(const real_t * x, int_t & cell_id)
{
  auto num_owned = num_owned_cells();

  if (num_dims_==1) {
    const auto & cell2verts = connectivity(1, 0);
    auto is_inside = find_cell_unstruct_1d(
        num_owned,
        cell2verts.offsets.data(),
        cell2verts.indices.data(),
        vertex_coords_.data(),
        x,
        cell_id);
    return is_inside;
  }
  else if (num_dims_==2) {
    const auto & cell2verts = connectivity(2, 0);
    auto is_inside = find_cell_unstruct_2d(
        num_owned,
        cell2verts.offsets.data(),
        cell2verts.indices.data(),
        vertex_coords_.data(),
        x,
        cell_id);
    return is_inside;
  }
  else {
    const auto & face2verts = connectivity(2, 0);
    const auto & cell2faces = connectivity(3, 2);
    auto is_inside = find_cell_unstruct_3d(
        num_owned,
        cell2faces.offsets.data(),
        cell2faces.indices.data(),
        face2verts.offsets.data(),
        face2verts.indices.data(),
        vertex_coords_.data(),
        x,
        cell_id);
    return is_inside;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Set num entities
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::set_num_vertices()
{
  auto nverts = vertex_coords_.size() / num_dims_;
  num_owned_entities(0) = nverts;
  num_all_entities(0) = nverts;
}

void mesh_unstruct_block_t::set_num_entities(int_t dim)
{
  auto it = connectivity_.find({dim, 0});
  if (it != connectivity_.end()) {
    auto n = it->second.size();
    num_owned_entities(dim) = n;
    num_all_entities(dim) = n;
  }
}
  
void mesh_unstruct_block_t::set_num_cells()
{
  auto ncells = local2global_.at(num_dims_).size();
  num_owned_entities(num_dims_) = ncells;
  num_all_entities(num_dims_) = ncells;
}
////////////////////////////////////////////////////////////////////////////////
/// Install boundaries
////////////////////////////////////////////////////////////////////////////////
bool mesh_unstruct_block_t::find_boundary_from_id(int_t bid, int_t & side_id) const
{
  auto bit = boundary_id2side_.find(bid);
  if (bit == boundary_id2side_.end()) return false;
  side_id = bit->second;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a boundary label
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::add_boundary_mapping(int_t side_id, int_t bid)
{
  boundary_side2id_[side_id].emplace_back(bid);
  boundary_id2side_.emplace(bid, side_id);
}
  
void mesh_unstruct_block_t::add_boundary_mappings(
    int_t side_id,
    const std::vector<int_t> labels)
{
  auto & blabels = boundary_side2id_[side_id];
  blabels.insert(blabels.end(), labels.begin(), labels.end());
  for (auto l : labels)
    boundary_id2side_.emplace(l, side_id);
}


////////////////////////////////////////////////////////////////////////////////
/// Install boundaries
////////////////////////////////////////////////////////////////////////////////
void mesh_unstruct_block_t::install_boundaries(
    const std::vector<std::unique_ptr<mesh_boundary_t>> & boundaries)
{
  boundary_side2id_.clear();
  boundary_id2side_.clear();

  for (const auto & bnd : boundaries)
    if (auto b = dynamic_cast<mesh_unstruct_boundary_t*>(bnd.get())) 
      add_boundary_mapping(b->side_id(), b->id());
}

  

} // namespace
