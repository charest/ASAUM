#include "mesh_struct_block.hpp"

#include "amr_flags.hpp"
#include "connect.hpp"
#include "ijk/block_ijk.hpp"
#include "mesh_boundary.hpp"
#include "mesh_struct_node.hpp"
#include "mesh_struct_block_cpu.hpp"
#include "comm/entity_info.hpp"
#include "math/subdivide.hpp"
#include "utils/cast.hpp"

#ifdef HAVE_EXODUS
#include "io/exo_writer.hpp"
#endif
#include "io/vtk_writer.hpp"

#include <cmath>

namespace prl {
  
////////////////////////////////////////////////////////////////////////////////
/// Helper functiion to remove unique entries
////////////////////////////////////////////////////////////////////////////////
template<typename T>
void remove_unique(T & vec)
{
  std::sort(vec.begin(), vec.end());
  auto last = std::unique(vec.begin(), vec.end());
  vec.erase(last, vec.end());
}

////////////////////////////////////////////////////////////////////////////////
/// Structured mesh block
////////////////////////////////////////////////////////////////////////////////
mesh_struct_block_t::mesh_struct_block_t(
    int_t region,
    const long_t * dims,
    int_t ndims) : 
  mesh_block_t(ndims),
  region_(region),
  dims_(dims, dims+ndims),
  iterator_(ndims, dims)
{

  // vertices
  auto vertices = iterator_.vertices();
  auto num_verts = vertices.size();

  vertex_coords_.clear();
  vertex_coords_.resize(num_dims_ * num_verts);

  num_owned_entities(0) = num_verts;
  num_all_entities(0) = num_verts;
  
  // faces
  int_t num_faces = 0;
  for (int_t d=0; d<num_dims_; ++d)
    num_faces += iterator_.faces(d).size();
  
  num_owned_entities(num_dims_-1) = num_faces;
  num_all_entities(num_dims_-1) = num_faces;
    
  // cells
  auto cells = iterator_.cells();
  auto num_cells = cells.size();

  num_owned_entities(num_dims_) = num_cells;
  num_all_entities(num_dims_) = num_cells;
  
}

////////////////////////////////////////////////////////////////////////////////
/// Structured mesh block
////////////////////////////////////////////////////////////////////////////////
mesh_struct_block_t::mesh_struct_block_t(const block_ijk_t & block) : 
  mesh_struct_block_t(block.label(), block.dims(), block.num_dims())
{

  auto block_coords = block.coords();

  auto vertices = iterator_.vertices();
  auto vit = vertices.begin();
  do {
    long_t is[] = {vit.i, vit.j, vit.k};
    block.spacing( &block_coords[0], is, &dims_[0], &vertex_coords_[vit.id()*num_dims_] );
  } while (vit.next());

}
  
////////////////////////////////////////////////////////////////////////////////
/// Cast to structured node
////////////////////////////////////////////////////////////////////////////////
mesh_struct_node_t * mesh_struct_block_t::my_node()
{ return dynamic_cast<mesh_struct_node_t *>(node_); }

const mesh_struct_node_t * mesh_struct_block_t::my_node() const
{ return dynamic_cast<const mesh_struct_node_t *>(node_); }


////////////////////////////////////////////////////////////////////////////////
/// Get neighbor info
////////////////////////////////////////////////////////////////////////////////
bool mesh_struct_block_t::has_neighbor_at(sector_t sec) const
{
  auto n = my_node();
  return n && n->neighbors().count(sec);
}

const std::vector<struct_neighbor_t> &
mesh_struct_block_t::neighbors_at(sector_t sec) const
{ return my_node()->neighbors().at(sec); }

const sector_map_t<std::vector<struct_neighbor_t>> &
mesh_struct_block_t::neighbors() const
{ return my_node()->neighbors(); }

const std::vector<struct_shared_t> &
mesh_struct_block_t::shared() const
{ return my_node()->shared(); }

////////////////////////////////////////////////////////////////////////////////
/// Count vertices i own
////////////////////////////////////////////////////////////////////////////////
int_t mesh_struct_block_t::count_vertices() const
{
  int_t my_verts = 0;
  
  const auto & permutations = sector_permutations[num_dims_-1];

  const auto & vertices = iterator_.vertices();
  auto vit = vertices.begin();

  auto my_id = id();
      
  do {
      
    auto block_owner = my_id;

    //--- boundary vertex
    if (vit.is_boundary()) {
      
      int_t pos[3];
      vit.pos(pos);
      
      for (const auto & sector : permutations) {
        auto is_sector_bnd = vit.is_sector_boundary(sector.begin());
        if (is_sector_bnd && has_neighbor_at(sector)) {
          auto neighs = neighbors_at(sector);
          for (const auto & n : neighs) {
            if (n.is_used())
              if (n.vert_is_inside(pos))
                block_owner = std::min(block_owner, n.id());
          }
        }
      }
    }
    
    if (block_owner == my_id) my_verts++;

  } while (vit.next());

  return my_verts;
}

////////////////////////////////////////////////////////////////////////////////
/// Number vertices
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::number_vertices(
    const std::vector<int_t> & block_offsets,
    long_t vert_start,
    std::vector<entity_info_t> & shared_info,
    std::vector<entity_info_t> & ghost_info)
{
  auto vertices = iterator_.vertices();
  auto num_verts = vertices.size();
  auto vit = vertices.begin();

  auto & vert_gids = global_ids(0);
  vert_gids.clear();
  vert_gids.resize(num_verts);

  auto vert_id = vert_start;
  std::vector<const struct_neighbor_t *> neighbors;

  const auto & permutations = sector_permutations[num_dims_-1];

  auto my_id = id();

  do {

    auto lid = vit.id();
      
    neighbors.clear();
    int_t block_owner = my_id;
    const struct_neighbor_t * neigh = nullptr;

    // boundary vertex
    if (vit.is_boundary()) {
      
      int_t pos[3];
      vit.pos(pos);
      
      for (const auto & sector : permutations) {

        auto is_sector_bnd = vit.is_sector_boundary(sector.begin());
        if (is_sector_bnd && has_neighbor_at(sector)) {
          const auto & neighs = neighbors_at(sector);
          for (const auto & n : neighs) {
            if (n.is_used()  && n.vert_is_inside(pos)) {
              auto nid = n.id();
              neighbors.emplace_back(&n);
              if (nid < block_owner) {
                block_owner = nid;
                neigh = &n;
              }
            } // inside

          } // neighs
        } // has_neighbor
      } // sector
    
    } // is_boundary

    // remove duplicates
    std::sort(neighbors.begin(), neighbors.end(),
        [](auto a, auto b) { return a->id()<b->id(); });

    auto last = std::unique(neighbors.begin(), neighbors.end(),
        [](auto & a, auto & b) { return a->id()==b->id(); });

    neighbors.erase(last, neighbors.end());

    //--- i own this vertex
    if (block_owner == my_id) {
      vert_gids[lid] = vert_id;
      vert_id++;
      for (auto n : neighbors) {
        auto nblk = n->id();
        auto nrnk = owner(block_offsets, nblk);
        shared_info.emplace_back(lid, nrnk, nblk, lid);
      }
    }
    //--- i need this vertex
    else {
      vert_gids[lid] = -1;
      int_t ipos[] = {vit.i, vit.j, vit.k};
      neigh->me2donor_verts(ipos);
      logical_block_t neigh_block(num_dims_, neigh->node()->as_struct()->dims());
      auto neigh_vertices = neigh_block.vertices();
      auto nid = neigh_vertices.id(ipos[0], ipos[1], ipos[2]);
      auto nrnk = owner(block_offsets, block_owner);
      ghost_info.emplace_back(lid, nrnk, block_owner, nid);
    }

  } while (vit.next());

}
  
////////////////////////////////////////////////////////////////////////////////
/// Add a boundary label
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::add_boundary_mapping(bool is_hi, int_t dir, int_t bid)
{
  boundary_side2id_[dir][is_hi].emplace_back(bid);
  boundary_id2side_.emplace(bid, label_pair_t{is_hi, dir});
}
  
void mesh_struct_block_t::add_boundary_mappings(
    bool is_hi,
    int_t dir,
    const std::vector<int_t> labels)
{
  auto & blabels = boundary_side2id_[dir][is_hi];
  blabels.insert(blabels.end(), labels.begin(), labels.end());
  for (auto l : labels)
    boundary_id2side_.emplace(l, label_pair_t{is_hi, dir});
}

////////////////////////////////////////////////////////////////////////////////
/// find a boundary label
////////////////////////////////////////////////////////////////////////////////
bool mesh_struct_block_t::find_boundary_from_id(int_t bid, bool & is_hi, int_t & dir) const
{
  auto bit = boundary_id2side_.find(bid);
  if (bit == boundary_id2side_.end()) return false;
  dir = bit->second.dim;
  is_hi = bit->second.is_hi;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
/// Install boundaries
////////////////////////////////////////////////////////////////////////////////
bool mesh_struct_block_t::install_boundaries(
    const block_ijk_t & block_ijk,
    const std::vector<std::unique_ptr<mesh_boundary_t>> & boundaries)
{
  for (int_t d=0; d<num_dims_; ++d)
    for (int_t s=0; s<2; ++s)
      boundary_side2id_[d][s].clear();
  boundary_id2side_.clear();
  
  for (const auto & bnd : boundaries) {
    if (auto b = dynamic_cast<mesh_struct_boundary_t*>(bnd.get())) {
      auto err = install_boundary(block_ijk, *b);
      if (err) {
        std::cout << "Boundary' " << bnd->name();
        std::cout << "' does not link to a face." << std::endl;
        return true;
      }
    }
  
  } // bounndaries

  return false;
}

////////////////////////////////////////////////////////////////////////////////
/// Add a boundary
////////////////////////////////////////////////////////////////////////////////
bool mesh_struct_block_t::install_boundary(
    const block_ijk_t & block_ijk,
    const mesh_struct_boundary_t & bnd)
{
  std::vector<const block_ijk_t*> leafs;

  // first get the root block
  auto rt = block_ijk.root();

  for (const auto & vs : bnd.patches()) {

    // see if this matches a sector in the root
    if (auto res = rt->sector(vs)) {

      const auto & sector = res.value;
      if (sector.nnz()!=1) return true;

      // now get the leaf blocks at this sector
      leafs.clear();
      rt->find_sectors(sector, leafs);
      
      std::sort(leafs.begin(), leafs.end());

      // Is this block part of this list, if so, it has this boundary
      if (std::binary_search(leafs.begin(), leafs.end(), &block_ijk)) {
        auto dim = sector.first_non_zero();
        auto dir = sector.is_hi(dim);
        add_boundary_mapping(dir, dim, bnd.id());
      }

    } // sector

  } // patches

  return false;
}
  
////////////////////////////////////////////////////////////////////////////////
/// build neighbors
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::build_neighbors()
{
  for (auto & neigh : desired_neighbors_)
    if (neighbors_.count(neigh) == 0)
      build_neighbors(neigh.first, neigh.second);
}

////////////////////////////////////////////////////////////////////////////////
/// compute neighbors
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::build_neighbors(int_t dim, int_t thru) 
{

  auto & neigh = neighbors_[{dim, thru}];
  auto & indices = neigh.indices;
  auto & offsets = neigh.offsets;

  indices.clear();
  offsets.clear();

  offsets.resize(num_all_cells()+1, 0);
        
  auto num_owned = num_owned_cells();

  //----------------------------------------------------------------------------
  // Cells
  if (dim == num_dims_) {
  
    //----------------------------------
    // Create a utilitiy function for
    // addig ghosts
      
    crs_t<int_t, int_t> ghost_ent_to_ghost_cell;

    auto ghost_cell_to_ghost_ent = ghost_connectivity_.find({num_dims_, thru});
    auto has_ghost = ghost_cell_to_ghost_ent != ghost_connectivity_.end();
    if (has_ghost)
      prl::transpose(
        ghost_cell_to_ghost_ent->second,
        ghost_ent_to_ghost_cell);

    auto cell_to_ghost_ent = cells2ghost_.find(thru);

    //--- The lambda
    auto add_ghost_cell = [&](const int_t id) {
      if (!has_ghost) return;

      const auto & c2ent = cell_to_ghost_ent->second;
      auto it = c2ent.find(id);
      if (it == c2ent.end()) return;
      
      auto before_ghost = indices.size();

      for (auto e : it->second)
        for (auto c :  ghost_ent_to_ghost_cell.at(e))
          indices.emplace_back(num_owned+c);

      auto first = indices.begin()+before_ghost;
      std::sort( first, indices.end() );
      auto last = std::unique( first, indices.end() );
      indices.erase(last, indices.end() );
    };
    
    //----------------------------------
    // Loop over cells

    auto cells = iterator_.cells();
    auto cit = cells.begin();
    do {
    
      auto im = cit.is_iboundary(false);
      auto ip = cit.is_iboundary(true);
      
      auto jm = cit.is_jboundary(false);
      auto jp = cit.is_jboundary(true);
      
      auto km = cit.is_kboundary(false);
      auto kp = cit.is_kboundary(true);
    
      //--- Cells -> Verts
      if (thru == 0) {
        if (num_dims_ == 1) {
          if (!im) indices.emplace_back( cit.offset_id(-1, 0, 0) );
          if (!ip) indices.emplace_back( cit.offset_id( 1, 0, 0) );
        }
        else if (num_dims_ == 2) {
          if (!im && !jm) indices.emplace_back( cit.offset_id(-1, -1, 0) );
          if (       !jm) indices.emplace_back( cit.offset_id( 0, -1, 0) );
          if (!ip && !jm) indices.emplace_back( cit.offset_id( 1, -1, 0) );
          if (!ip       ) indices.emplace_back( cit.offset_id( 1,  0, 0) );
          if (!ip && !jp) indices.emplace_back( cit.offset_id( 1,  1, 0) );
          if (       !jp) indices.emplace_back( cit.offset_id( 0,  1, 0) );
          if (!im && !jp) indices.emplace_back( cit.offset_id(-1,  1, 0) );
          if (!im       ) indices.emplace_back( cit.offset_id(-1,  0, 0) );
        }
        else {
          if (              !km) indices.emplace_back( cit.offset_id( 0,  0, -1) );
          if (!im && !jm && !km) indices.emplace_back( cit.offset_id(-1, -1, -1) );
          if (       !jm && !km) indices.emplace_back( cit.offset_id( 0, -1, -1) );
          if (!ip && !jm && !km) indices.emplace_back( cit.offset_id( 1, -1, -1) );
          if (!ip &&        !km) indices.emplace_back( cit.offset_id( 1,  0, -1) );
          if (!ip && !jp && !km) indices.emplace_back( cit.offset_id( 1,  1, -1) );
          if (       !jp && !km) indices.emplace_back( cit.offset_id( 0,  1, -1) );
          if (!im && !jp && !km) indices.emplace_back( cit.offset_id(-1,  1, -1) );
          if (!im &&        !km) indices.emplace_back( cit.offset_id(-1,  0, -1) );
          
          //if (                 ) indices.emplace_back( cit.offset_id( 0,  0,  0) );
          if (!im && !jm       ) indices.emplace_back( cit.offset_id(-1, -1,  0) );
          if (       !jm       ) indices.emplace_back( cit.offset_id( 0, -1,  0) );
          if (!ip && !jm       ) indices.emplace_back( cit.offset_id( 1, -1,  0) );
          if (!ip              ) indices.emplace_back( cit.offset_id( 1,  0,  0) );
          if (!ip && !jp       ) indices.emplace_back( cit.offset_id( 1,  1,  0) );
          if (       !jp       ) indices.emplace_back( cit.offset_id( 0,  1,  0) );
          if (!im && !jp       ) indices.emplace_back( cit.offset_id(-1,  1,  0) );
          if (!im              ) indices.emplace_back( cit.offset_id(-1,  0,  0) );
          
          if (              !kp) indices.emplace_back( cit.offset_id( 0,  0,  1) );
          if (!im && !jm && !kp) indices.emplace_back( cit.offset_id(-1, -1,  1) );
          if (       !jm && !kp) indices.emplace_back( cit.offset_id( 0, -1,  1) );
          if (!ip && !jm && !kp) indices.emplace_back( cit.offset_id( 1, -1,  1) );
          if (!ip &&        !kp) indices.emplace_back( cit.offset_id( 1,  0,  1) );
          if (!ip && !jp && !kp) indices.emplace_back( cit.offset_id( 1,  1,  1) );
          if (       !jp && !kp) indices.emplace_back( cit.offset_id( 0,  1,  1) );
          if (!im && !jp && !kp) indices.emplace_back( cit.offset_id(-1,  1,  1) );
          if (!im &&        !kp) indices.emplace_back( cit.offset_id(-1,  0,  1) );
        }

      }
      
      //--- Cells -> Edges
      else if (thru == 1) {

        if (num_dims_ == 2) {
          if (!im) indices.emplace_back( cit.offset_id( -1,  0,  0) );
          if (!jm) indices.emplace_back( cit.offset_id(  0, -1,  0) );
          if (!ip) indices.emplace_back( cit.offset_id(  1,  0,  0) );
          if (!jp) indices.emplace_back( cit.offset_id(  0,  1,  0) );
        }
        else if (num_dims_ == 3) {
          if (              !km) indices.emplace_back( cit.offset_id( 0,  0, -1) );
          if (       !jm && !km) indices.emplace_back( cit.offset_id( 0, -1, -1) );
          if (!ip &&        !km) indices.emplace_back( cit.offset_id( 1,  0, -1) );
          if (       !jp && !km) indices.emplace_back( cit.offset_id( 0,  1, -1) );
          if (!im &&        !km) indices.emplace_back( cit.offset_id(-1,  0, -1) );
          
          //indices.emplace_back( cit.offset_id( 0,  0,  0) );
          if (!im && !jm       ) indices.emplace_back( cit.offset_id(-1, -1,  0) );
          if (       !jm       ) indices.emplace_back( cit.offset_id( 0, -1,  0) );
          if (!ip && !jm       ) indices.emplace_back( cit.offset_id( 1, -1,  0) );
          if (!ip              ) indices.emplace_back( cit.offset_id( 1,  0,  0) );
          if (!ip && !jp       ) indices.emplace_back( cit.offset_id( 1,  1,  0) );
          if (       !jp       ) indices.emplace_back( cit.offset_id( 0,  1,  0) );
          if (!im && !jp       ) indices.emplace_back( cit.offset_id(-1,  1,  0) );
          if (!im              ) indices.emplace_back( cit.offset_id(-1,  0,  0) );
          
          if (              !kp) indices.emplace_back( cit.offset_id( 0,  0,  1) );
          if (       !jm && !kp) indices.emplace_back( cit.offset_id( 0, -1,  1) );
          if (!ip &&        !kp) indices.emplace_back( cit.offset_id( 1,  0,  1) );
          if (       !jp && !kp) indices.emplace_back( cit.offset_id( 0,  1,  1) );
          if (!im &&        !kp) indices.emplace_back( cit.offset_id(-1,  0,  1) );
        }
      }
      //--- Cells -> Faces
      else if (thru == 2) {
          
        if (!im) indices.emplace_back( cit.offset_id( -1,  0,  0) );
        if (!ip) indices.emplace_back( cit.offset_id(  1,  0,  0) );
        if (!jm) indices.emplace_back( cit.offset_id(  0, -1,  0) );
        if (!jp) indices.emplace_back( cit.offset_id(  0,  1,  0) );
        if (!km) indices.emplace_back( cit.offset_id(  0,  0, -1) );
        if (!kp) indices.emplace_back( cit.offset_id(  0,  0,  1) );

      }

      //--- Check ghost
      auto cid = cit.id();
      add_ghost_cell(cid);

      offsets[cid+1] = indices.size();

    } while(cit.next());

  }
  
  //----------------------------------------------------------------------------
  // Edges
  else if (dim==1) {
    THROW_ERROR("HAVE NOTT IMPLEMENTTEED EDGE NEIGHBORS");
  }
  
  //----------------------------------------------------------------------------
  // Faces
  else if (dim==2) {
    THROW_ERROR("HAVE NOTT IMPLEMENTTEED FACE NEIGHBORS");
  }
  

}

////////////////////////////////////////////////////////////////////////////////
/// build connecitivy
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::build_connectivity()
{
  for (auto & conn : desired_connectivity_)
    if (connectivity_.count(conn) == 0)
      build_connectivity(conn.first, conn.second);
}


////////////////////////////////////////////////////////////////////////////////
/// compute face neighbors
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::build_connectivity(int_t a, int_t b) 
{

  auto max = std::max(a, b);
  auto min = std::min(a, b);

  //----------------------------------------------------------------------------
  // Cells to faces
  //----------------------------------------------------------------------------
  if (max == num_dims_ && min == num_dims_-1) {
  
    logical_face_t faces[3] = {
      iterator_.faces(0),
      iterator_.faces(1),
      iterator_.faces(2)};
    auto cells = iterator_.cells();
    auto cit = cells.begin();

    crs_t<int_t, int_t> conn;
    auto & indices = conn.indices;
    auto & offsets = conn.offsets;

    auto num_cells = num_all_cells();
    offsets.reserve( num_cells+1 );
    indices.reserve( num_cells*std::pow(2, num_dims_) );
    
    offsets.emplace_back(0);

    //--- internal cells
    do {
     
      for (int_t dir=0; dir<num_dims_; ++dir) {
        for (int_t side=0, off=-1; side<2; ++side, off+=2) {    
          auto f = cit.face_id(faces[dir], side);
          indices.emplace_back(f);
        }
      }

      offsets.emplace_back(indices.size());
    
    } while (cit.next());

    //--- ghost cells
    auto it = ghost_connectivity_.find({max, min});
    if (it != ghost_connectivity_.end()) {

      const auto & ghost2ent =  ghost2entities_.at(min);

      for (const auto & ents : it->second) {
        indices.reserve( indices.size() + ents.size() );
        for (auto e : ents) {
          auto res = ghost2ent.find(e);
          if (res != ghost2ent.end())
            indices.emplace_back(res->second);
        }
        offsets.emplace_back(indices.size());
      }
    }

    //--- transpose if needed
    if (a != max)
      prl::transpose(conn, connectivity_[{a,b}]);
    else
      connectivity_[{a,b}] = std::move(conn);

  }
  //----------------------------------------------------------------------------
  // Cells to vertex
  //----------------------------------------------------------------------------
  else if( max==num_dims_ && min== 0 )
  {
    
    auto nx = dims_[0];
    auto ny = dims_[1];
    auto nz = dims_[2];
    
    const auto & vertices = iterator_.vertices();
    const auto & cells = iterator_.cells();
    
    int num_cells = cells.size();
    int size_offset=pow(2,num_dims_);
    std::vector<int_t> cell2verts(size_offset * num_cells);

    //------------------------------------
    if (num_dims_ == 1) {
      
      for (int_t i=0; i<nx; ++i) {
        cell2verts[ i*2    ] = i;//left
        cell2verts[ i*2 + 1] = i+1;//right
      }
      
    }
    //------------------------------------
    else if (num_dims_ == 2) {
      
      for (int_t j=0; j<ny; ++j) {
        for (int_t i=0; i<nx; ++i) {
          auto aa = vertices.id(i+0, j+0, 0);
          auto bb = vertices.id(i+1, j+0, 0);
          auto cc = vertices.id(i+1, j+1, 0);
          auto dd = vertices.id(i+0, j+1, 0);
          auto id = cells.id(i, j, 0);

          cell2verts[ id*4    ] = aa;//bottom-left
          cell2verts[ id*4 + 1] = bb;//bottom-right
          cell2verts[ id*4 + 2] = cc;//top-right
          cell2verts[ id*4 + 3] = dd;//top-left
      }
    }

    }
    //------------------------------------
    else {
      for (int_t k=0; k<nz; ++k) {
        for (int_t j=0; j<ny; ++j) {
          for (int_t i=0; i<nx; ++i) {
            auto a = vertices.id(i+0, j+0, k+0);//-x-y-z
            auto b = vertices.id(i+1, j+0, k+0);//+x-y-z
            auto c = vertices.id(i+1, j+1, k+0);//+x+y-z
            auto d = vertices.id(i+0, j+1, k+0);//-x+y-z

            auto e = vertices.id(i+0, j+0, k+1);//-x-y+z
            auto f = vertices.id(i+1, j+0, k+1);//+x-y-z
            auto g = vertices.id(i+1, j+1, k+1);//+x+y+z
            auto h = vertices.id(i+0, j+1, k+1);//-x+y-z
            auto id = cells.id(i, j, k);
            cell2verts[ id*8    ] = a;
            cell2verts[ id*8 + 1] = b;
            cell2verts[ id*8 + 2] = c;
            cell2verts[ id*8 + 3] = d;
            cell2verts[ id*8 + 4] = e;
            cell2verts[ id*8 + 5] = f;
            cell2verts[ id*8 + 6] = g;
            cell2verts[ id*8 + 7] = h;
          }
        }
      }
    
    }
    // now fill in.
    crs_t<int_t, int_t> conn;
    auto & indices = conn.indices;
    auto & offsets = conn.offsets;
    indices.resize(num_cells);
    indices.assign( cell2verts.begin(), cell2verts.end());
    offsets.resize(num_cells+1, 0);
    
    for(int i =0;i<num_cells;++i)
    {
      offsets[i+1] = offsets[i]+ size_offset;
    }//i
    
    //--- transpose if needed
    if (a != max)
      prl::transpose(conn, connectivity_[{a,b}]);
    else
      connectivity_[{a,b}] = std::move(conn);

  }
  //----------------------------------------------------------------------------
  // Unsupported
  //----------------------------------------------------------------------------
  else {

    THROW_ERROR("Connectivity from " << a << " to " << b << " for structured"
        << " meshes is not implemented yet.");

  }
}

////////////////////////////////////////////////////////////////////////////////
/// builld geeometry
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::_build_geometry(bool with_face_geom, bool with_cell_geom)
{
  if (with_face_geom)
    struct_face_geometry(
        iterator_,
        vertex_coords_.data(),
        face_centroids_.data(),
        face_normals_.data(),
        face_areas_.data());

  if (with_cell_geom) {
    struct_cell_geometry(
        iterator_,
        vertex_coords_.data(),
        cell_centroids_.data(),
        cell_volumes_.data());
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Extract the mesh boundary
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::extract_boundary(
    int_t boundary_id,
    std::vector<boundary_face_t> & bfaces)
{
  int_t ibeg[3] = {0, 0, 0};
  int_t iend[3] = {1, 1, 1};
  int_t ioff[3] = {0, 0, 0};

  const auto & cells = iterator_.cells();

  
  // check all boundaries
  bool side;
  int_t dim;
  auto has_boundary = find_boundary_from_id(boundary_id, side, dim);
  if (!has_boundary) return;

  int_t dir = side ? 1:-1;
        
  //--- set offsets
  for (int_t i=0; i<num_dims_; ++i) {
    ibeg[i] = 0;
    iend[i] = dims_[i];
    ioff[i] = 0;
  }
  ibeg[dim] = side==0 ? 0 : dims_[dim];
  iend[dim] = ibeg[dim]+1;
  ioff[dim] = side==0 ? 0 : -1;
  
  const auto & faces = iterator_.faces(dim);
  
  //--- loop over boundaries
  for (int_t k=ibeg[2]; k!=iend[2]; ++k)
    for (int_t j=ibeg[1]; j<iend[1]; ++j)
      for (int_t i=ibeg[0]; i<iend[0]; ++i) {
        auto fid = faces.id(i, j, k);
        auto cid = cells.id(i+ioff[0], j+ioff[1], k+ioff[2]);
        bfaces.emplace_back(fid, dir, cid);
      }
  
}


////////////////////////////////////////////////////////////////////////////////
/// vtk output
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::output(vtk_writer_t & vtk)  const
{
  vtk.write_points(vertex_coords_, num_dims_);
}

////////////////////////////////////////////////////////////////////////////////
/// exo output
////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_EXODUS
void mesh_struct_block_t::output(
    exo_writer_t & exo,
    const std::string & label,
    const std::vector<side_set_info_t> &)  const
{
  auto num_cells = num_owned_cells();
  if (num_cells==0) return;

  auto num_verts = num_owned_vertices();

  // write param:
  auto params = exo_writer_t::make_params(label.c_str());
  params.num_dim = num_dims_;
  params.num_nodes = num_verts;
  params.num_elem_blk = 1; // only one region for now
  params.num_elem = num_cells; 
  exo.write_params(params);
  
  // check the integer type used in the exodus file
  auto int64 = exo.is_int64();
  
  //----------------------------------------------------------------------------
  // Determine cell connectivity
      
  auto verts_per_cell = std::pow(2, num_dims_);
  
  std::vector<int_t> cell2verts(verts_per_cell * num_cells);
  
  shape_t cell_type;

  auto nx = dims_[0];
  auto ny = dims_[1];
  auto nz = dims_[2];

  const auto & vertices = iterator_.vertices();
  const auto & cells = iterator_.cells();
  
  //------------------------------------
  if (num_dims_ == 1) {

    cell_type = shape_t::line;

    for (int_t i=0; i<nx; ++i) {
      cell2verts[ i*2    ] = i;
      cell2verts[ i*2 + 1] = i+1;
    }

  }
  //------------------------------------
  else if (num_dims_ == 2) {
    
    cell_type = shape_t::quad;

    for (int_t j=0; j<ny; ++j) {
      for (int_t i=0; i<nx; ++i) {
        auto a = vertices.id(i+0, j+0, 0);
        auto b = vertices.id(i+0, j+1, 0);
        auto c = vertices.id(i+1, j+1, 0);
        auto d = vertices.id(i+1, j+0, 0);
        auto id = cells.id(i, j, 0);
        cell2verts[ id*4    ] = a;
        cell2verts[ id*4 + 1] = b;
        cell2verts[ id*4 + 2] = c;
        cell2verts[ id*4 + 3] = d;
      }
    }

  }
  //------------------------------------
  else {
    
    cell_type = shape_t::hex;

    for (int_t k=0; k<nz; ++k) {
      for (int_t j=0; j<ny; ++j) {
        for (int_t i=0; i<nx; ++i) {
          auto a = vertices.id(i+0, j+0, k+0);
          auto b = vertices.id(i+0, j+1, k+0);
          auto c = vertices.id(i+1, j+1, k+0);
          auto d = vertices.id(i+1, j+0, k+0);
          auto e = vertices.id(i+0, j+0, k+1);
          auto f = vertices.id(i+0, j+1, k+1);
          auto g = vertices.id(i+1, j+1, k+1);
          auto h = vertices.id(i+1, j+0, k+1);
          auto id = cells.id(i, j, k);
          cell2verts[ id*8    ] = a;
          cell2verts[ id*8 + 1] = b;
          cell2verts[ id*8 + 2] = c;
          cell2verts[ id*8 + 3] = d;
          cell2verts[ id*8 + 4] = e;
          cell2verts[ id*8 + 5] = f;
          cell2verts[ id*8 + 6] = g;
          cell2verts[ id*8 + 7] = h;
        }
      }
    }

  } // num_dims
  
  //----------------------------------------------------------------------------
  // Write cells

  // helper function for getting cell vertices
  auto cell_vertices_func = [&](auto c, auto & vert_list) {
    auto start = &cell2verts[ c * verts_per_cell ];
    auto end = &cell2verts[ (c+1) * verts_per_cell ];
    vert_list.insert(vert_list.end(), start, end);
  };
  
  // add the whole element block
  std::stringstream ss;
  ss << "Block " << id();
  
  auto type_str = exo_writer_t::geom_type(cell_type);
  auto region_id = 1;

  if (int64) {
    exo.write_element_block<long long>(
        region_id,
        ss.str().c_str(),
        type_str,
        num_cells,
        cell_vertices_func);
  }
  else {
    exo.write_element_block<int>(
        region_id,
        ss.str().c_str(),
        type_str,
        num_cells,
        cell_vertices_func);
  }
  
  //----------------------------------------------------------------------------
  // write final element mapping


  const auto & cell_global_ids = global_ids(num_dims_);
  auto cell_local2global_ref = make_array_ref(cell_global_ids.data(), num_cells);
  
  if (int64)
    exo.write_element_map<long long>(cell_local2global_ref);
  else
    exo.write_element_map<int>(cell_local2global_ref);
      
  //----------------------------------------------------------------------------
  // write coordinates
  std::vector<real_t> coordinates(num_dims_ * num_verts);

  for (int_t v=0; v<num_verts; ++v) {
    for (int_t d=0; d<num_dims_; ++d )
      coordinates[v + d*num_verts] = vertex_coords_[v*num_dims_ + d];
  }
  exo.write_point_coords(coordinates);
  
  //----------------------------------------------------------------------------
  // write final vertex mapping
  
  const auto & vertex_global_ids = global_ids(0);
  auto vertex_local2global_ref = make_array_ref(vertex_global_ids.data(), num_verts);

  if (int64)
    exo.write_node_map<long long>(vertex_local2global_ref);
  else
    exo.write_node_map<int>(vertex_local2global_ref);
} // output
#endif

////////////////////////////////////////////////////////////////////////////////
/// find a cell
////////////////////////////////////////////////////////////////////////////////
bool mesh_struct_block_t::find_cell(const real_t * x, int_t & cell_id)
{
  auto is_inside = find_cell_struct(
      iterator_,
      vertex_coords_.data(),
      x,
      cell_id);
  return is_inside;
}

////////////////////////////////////////////////////////////////////////////////
/// pack a cell into a buffer
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::pack(
    int_t local_id,
    std::vector<byte_t>& buffer,
    bool /*as_ghost*/) const
{
  
  //--------------------------------------------------------------------------
  // Determine the cell
  const auto & cells = iterator_.cells();
  auto cit = cells.at(local_id);
  
  //--------------------------------------------------------------------------
  // pack cell info
  
  // add mesh global id to buffer
  long_t global_id = global_ids(num_dims_)[local_id];
  cast_insert( &global_id, 1, buffer );
  
  //--------------------------------------------------------------------------
  // pack vertices
  
  int_t vs[8];
  int_t num_verts = std::pow(2, num_dims_);
  const auto & vertices = iterator_.vertices();
  cit.vertex_ids(vertices, vs);

  const auto & vert_gids = global_ids(0);
    
  // add num verts to buffer
  cast_insert( &num_verts, 1, buffer );
  // add global vertex indices to buffer
  for (int_t v=0; v<num_verts; ++v) {
    long_t gid = vert_gids[vs[v]];
    cast_insert( &gid, 1, buffer );
  }

  //--------------------------------------------------------------------------
  // pack faces
  if (num_dims_ > 2) {
  
    int_t num_faces = 2*num_dims_;
    cast_insert( &num_faces, 1, buffer );
    
    num_verts = std::pow(2, num_dims_-1);

    for (int_t dim=0; dim<num_dims_; ++dim) {
      const auto & face = iterator_.faces(dim);
      for (int_t side=0; side<2; ++side) {

        auto fit = cit.face(face, side);
        fit.vertex_ids(vertices, vs);

        // add num verts to buffer
        cast_insert( &num_verts, 1, buffer );

        // add global vertex indices to buffer
        for (int_t v=0; v<num_verts; ++v) {
          long_t gid = vert_gids[vs[v]];
          cast_insert( &gid, 1, buffer );
        }

      } // side
    } // dim

  } // 3d
  
}

////////////////////////////////////////////////////////////////////////////////
/// unpack a cell
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::unpack(const byte_t * & buffer, bool /*as_ghost*/)
{
  
  //--------------------------------------------------------------------------
  // Cell
  
  // get global mesh id
  long_t global_id;
  uncast( buffer, 1, &global_id );
  
  // add mapping
  add_ghost_entity(global_id, num_dims_);

  //--------------------------------------------------------------------------
  // Vertices
  
  // get the number of vertices
  int_t num_verts;
  uncast( buffer, 1, &num_verts );
    
  std::vector<long_t> vs(num_verts);
  uncast(buffer, num_verts, vs.data());

  std::vector<int_t> new_vs( num_verts );
  for (int_t v=0; v<num_verts; ++v)
    new_vs[v] = add_ghost_entity(vs[v], 0);

  // add vertex connecitivity
  auto & cell2vertices = ghost_connectivity_[{num_dims_, 0}]; 
  cell2vertices.push_back(new_vs.begin(), new_vs.end());
  
  //--------------------------------------------------------------------------
  // 2D Edges
  if (num_dims_ == 2) {
    
    std::vector<int_t> es(num_verts);

    int_t edge_vs[2];
    for (int_t v=0; v<num_verts; ++v) {
      edge_vs[0] = new_vs[v];
      edge_vs[1] = new_vs[v+1==num_verts ? 0 : v+1];
      es[v] = add_ghost_entity(&edge_vs[0], &edge_vs[2], 1);
    }

    auto & cell2edges = ghost_connectivity_[{2, 1}];
    cell2edges.push_back(es.begin(), es.end());

  }
  
  //--------------------------------------------------------------------------
  // 3D Faces
  else if (num_dims_ == 3) {

    const auto & vert_global2local = ghost_global2local_[0];
    
    auto & cell2edges = ghost_connectivity_[{3, 1}];
    auto & cell2faces = ghost_connectivity_[{3, 2}];

    // get the number of faces
    int_t num_faces;
    uncast( buffer, 1, &num_faces );

    std::vector<int_t> fs;
    fs.reserve(num_faces);
    
    std::vector<int_t> es;
    int_t edge_vs[2];

    for (int_t f=0; f<num_faces; ++f) {

      // get the number of vertices
      int_t num_verts;
      uncast( buffer, 1, &num_verts );
      
      // retrieve global vertex ids
      vs.clear();
      vs.resize(num_verts);
      uncast( buffer, vs.size(), vs.data() );

      // convert vertex global ids to local ids
      new_vs.clear();
      new_vs.reserve(num_verts);

      for (auto & v : vs) {
        auto it = vert_global2local.find(v);
        assert(it != vert_global2local.end() && "could not find vertex");
        new_vs.emplace_back(it->second);
      }

      // add facce
      fs.emplace_back( add_ghost_entity(&new_vs[0], &new_vs[num_verts], 2) );
      
      // add edges
      for (int_t v=0; v<num_verts; ++v) {
        edge_vs[0] = new_vs[v];
        edge_vs[1] = new_vs[v+1==num_verts ? 0 : v+1];
        es.emplace_back( add_ghost_entity(&edge_vs[0], &edge_vs[1], 1) );
      }
  
    } 
    
    remove_unique(es);

    // now append to the connectivity list
    cell2faces.push_back( fs.begin(), fs.end() );
    cell2edges.push_back( es.begin(), es.end() );
      

  } // 3D

}

////////////////////////////////////////////////////////////////////////////////
/// Add a ghost vertex
////////////////////////////////////////////////////////////////////////////////
int_t mesh_struct_block_t::add_ghost_entity(long_t gid, int_t dim)
{
  auto & local2global = ghost_local2global_[dim];
  auto & global2local = ghost_global2local_[dim];

  auto res = global2local.emplace(gid, global2local.size());
  auto id = res.first->second;
  if (res.second) local2global.emplace_back(gid);
  return id;
}
    
////////////////////////////////////////////////////////////////////////////////
/// Add a ghost face/edge
////////////////////////////////////////////////////////////////////////////////
int_t mesh_struct_block_t::add_ghost_entity(const int_t * begin, const int_t * end, int_t dim)
{
  auto & sorted_vertices2face = sorted_verts2entities_[dim];
  auto & face2vertices = ghost_connectivity_[{dim, 0}]; 

  // sort vertices for search
  std::vector<int_t> sorted_vs(begin, end);
  std::sort(sorted_vs.begin(), sorted_vs.end());
    
  // search for the face to see if it exists
  auto res = sorted_vertices2face.emplace( sorted_vs, sorted_vertices2face.size() );

  // no match, add it to the connectivity list
  if (!res.second) face2vertices.push_back( begin, end );

  return res.first->second;
}

////////////////////////////////////////////////////////////////////////////////
/// Add an intreblock face
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::add_intrablock_face(
    const face_iterator_t & f,
    const int_t dir,
    const int_t owner_id,
    const int_t ghost_id)
{
  auto fid = f.id();
  intrablock_faces_.emplace_back(fid, dir, owner_id, ghost_id);

  //------------------------------------
  // Vertices
  
  const int_t nverts_per_face = std::pow(2, num_dims_-1);
  
  // get face vertices
  int_t vs[4];
  const auto & vertices = iterator_.vertices();
  f.vertex_ids(vertices, vs);

  // use the vertex global id to assign a new unique id
  const auto & vert_gids = global_ids(0);
  auto & verts2ghost = entities2ghost_[0];
  auto & ghost2verts = ghost2entities_[0];
  for (int_t v=0; v<nverts_per_face; ++v) {
    auto lid = vs[v];
    auto gid = vert_gids[lid];
    vs[v] = add_ghost_entity(gid, 0);
    verts2ghost.emplace(lid, vs[v]);
    ghost2verts.emplace(vs[v], lid);
  }

  // add link to intrablock face
  auto & cells2verts = cells2ghost_[0];
  auto & these_verts = cells2verts[owner_id];
  these_verts.insert(these_verts.end(), &vs[0], &vs[nverts_per_face]);
  remove_unique(these_verts);
  
  //------------------------------------
  // 2D edge
  if (num_dims_==2) {
    
    //--- Add edge
    int_t edge_vs[] = { vs[0], vs[1] };
    auto eid_ghost = add_ghost_entity(&edge_vs[0], &edge_vs[2], 1);
  
    entities2ghost_[1].emplace(fid, eid_ghost);
    ghost2entities_[1].emplace(eid_ghost, fid);
      
    // add link to intrablock edge
    auto & cells2edges = cells2ghost_[1];
    auto & these_edges = cells2edges[owner_id];
    these_edges.push_back(eid_ghost);
    remove_unique(these_edges);

  }

  //------------------------------------
  // 3D face
  else if (num_dims_==3) {
    
    //--- Add edges

    int_t es[4];
    auto & edges2ghost = entities2ghost_[1];
    auto & ghost2edges = ghost2entities_[1];

    for (int_t dim=0, e=0; dim<num_dims_; ++dim) {
      if (dim == f.dim()) continue;

      edge_iterator_t edge_its[2];
      
      auto edges = iterator_.edges(dim);
      f.edges(edges, edge_its);

      // two edges per dimesion
      for (int_t i=0; i<2; ++i) {
        auto eid = edge_its[i].id();
      
        int_t edge_vs[2];
        edge_its[i].vertex_ids(vertices, edge_vs); 
        
        for (auto & v : edge_vs) v = verts2ghost.at(v);
        
        es[e] = add_ghost_entity(&edge_vs[0], &edge_vs[2], 1);

        edges2ghost.emplace(eid, es[e]);
        ghost2edges.emplace(es[e], eid);
        
        e++;
      }

    }
      
    // add link to intrablock edge
    auto & cells2edges = cells2ghost_[1];
    auto & these_edges = cells2edges[owner_id];
    these_edges.insert(these_edges.end(), es, es+nverts_per_face);
    remove_unique(these_edges);
    
    //--- Add face
    auto fid_ghost = add_ghost_entity(vs, vs+nverts_per_face, 2);
    entities2ghost_[2].emplace(fid, fid_ghost);
    ghost2entities_[2].emplace(fid_ghost, fid);

    // add link to intrablock face
    auto & cells2faces = cells2ghost_[2];
    auto & these_faces = cells2faces[owner_id];
    these_faces.push_back(fid_ghost);
    remove_unique(these_faces);
    
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Build a halo
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::build_halo(
    int_t num_ghost,
    bool with_corners,
    const std::vector<int_t> & block_offsets,
    const std::vector<long_t> & cell_dist,
    std::vector<entity_info_t> & shared_info,
    std::vector<entity_info_t> & ghost_info)
{
  // alias the data and clear it
  auto & gids = global_ids(num_dims_);

  // make sure there are at least as many cells in each direction as ghost
  // cells
  for (int_t d=0; d<num_dims_; ++d) {
    auto n = dims_[d];
    if (n < num_ghost) {
      THROW_ERROR("Asking for a halo of depth of " << num_ghost
          << ",  but block with id " << id() << " only has "
          << n << " cells in the " << d << "th dimension.");
    }
  }
  
  //----------------------------------------------------------------------------
  // count how many ghost cells
  int_t ncells_ghost = 0;

  for (const auto & neigh_pair : neighbors()){
    const auto & sector =  neigh_pair.first;
    const auto & sector_neighbors = neigh_pair.second;
    auto sector_dims = sector.nnz();

    if ( sector_neighbors.size() &&
        (sector_dims == 1 || with_corners) )
    {
      int_t nghost = 1;
      for (int_t d=0; d<num_dims_; ++d) {
        if (sector[d]==0)
          nghost *= dims_[d];
        else
          nghost *= num_ghost;
      }

      ncells_ghost += nghost;
    }

  } // sector

  // create space
  const auto & cells = iterator_.cells();
  auto ncells_owned = cells.size();
  auto ncells_all = ncells_owned + ncells_ghost;

  num_all_entities(num_dims_) = ncells_all;

  gids.resize(ncells_all);
  ghost_info.reserve(ncells_ghost);
  shared_info.reserve(ncells_ghost);

  // reset ghost
  intrablock_faces_.clear();
  intrablock_cell_to_cell_.clear();
  intrablock_cell_to_intrablock_faces_.clear();
  sorted_verts2entities_.clear();
  cells2ghost_.clear();
  entities2ghost_.clear();
  ghost2entities_.clear();
  ghost_connectivity_.clear();
  ghost_local2global_.clear();
  ghost_global2local_.clear();

  intrablock_faces_.reserve(ncells_ghost);
  
  //----------------------------------------------------------------------------
  // add ghost info
  
  // ghost start after owned
  auto cell_cnt = ncells_owned;
  
  // compute some dimmensions
  logical_face_t faces[] = {
     iterator_.faces(0),
     iterator_.faces(1),
     iterator_.faces(2)
  };

  for (const auto & neigh_pair : neighbors()) {
    const auto & sector =  neigh_pair.first;
    const auto & sector_neighbors = neigh_pair.second;
    auto sector_dims = sector.nnz();

    if (!with_corners && sector_dims>1) continue;

    auto first_nz = sector.first_non_zero();
    const auto & my_range = my_node()->range(sector);
            
    //--- Transform coordinates
    for (const auto & neigh : sector_neighbors) {

      if (!neigh.is_used()) continue;
    
      // compute ghost ranges
      const auto & neigh_range = neigh.range();
      auto gbeg = neigh_range.begin;
      auto gend = neigh_range.end;
      int_t goff[3] = {1, 1, 1};
    
      // Compute shared ranges
      auto sbeg = neigh_range.begin;
      int_t soff[3] = {1, 1, 1};
      
      for (int_t d=0; d<num_dims_; ++d) {
        if (sector[d] != 0) {
          gbeg[d] += sector[d];
          gend[d] = gbeg[d] + sector[d]*num_ghost;
          goff[d] = sector[d];
          soff[d] = -sector[d];
          sbeg[d] = my_range.begin[d];
        }
      }

      // rank and block id
      auto nnode = neigh.node()->as_struct();
      auto nblk = neigh.id();
      auto nrnk = owner(block_offsets, nblk);

      //transformation matrix and iterator
      auto me2neigh = neigh.tm_me2donor_cells();

      // neigbbor root block
      auto ncells = logical_block_t(num_dims_, nnode->dims()).cells();

      // detrmine neighbor id
      for (int_t kg=gbeg[2], ks=sbeg[2]; kg!=gend[2]; kg+=goff[2], ks+=soff[2]) {
        for (int_t jg=gbeg[1], js=sbeg[1]; jg!=gend[1]; jg+=goff[1], js+=soff[1]) {
          for (int_t ig=gbeg[0], is=sbeg[0]; ig!=gend[0]; ig+=goff[0], is+=soff[0]) {
            // transform ghost point
            int_t ipos[] = {ig, jg, kg};
            me2neigh.transform(ipos);
            // add ghost
            auto nid = ncells.id(ipos[0], ipos[1], ipos[2]);
            ghost_info.emplace_back(cell_cnt, nrnk, nblk, nid);
            // determine global id
            gids[cell_cnt] = cell_dist[nblk] + nid;
            // add boundary face (only on the first layer of ghost)
            int_t face_pos[] = {is, js, ks};
            if (sector_dims == 1 && face_pos[first_nz]==sbeg[first_nz]) {
              if (sector[first_nz]>0) face_pos[first_nz]++;
              const auto & face = faces[first_nz];
              auto f = face.at(face_pos[0], face_pos[1], face_pos[2]);
              auto lid = cells.id(is, js, ks);
              add_intrablock_face(f, sector[first_nz], lid, cell_cnt);
            }
            // increment
            cell_cnt++;
          }
        }
      }

    } // neighbors
  } // sector neighbors
  
  //----------------------------------------------------------------------------
  // add shared info
  
  for (const auto & neigh : shared()) {
    const auto & sector =  neigh.sector();
    auto sector_dims = sector.nnz();

    if (!neigh.is_used()) continue;
    if (!with_corners && sector_dims>1) continue;

    // Compute shared ranges
    const auto & neigh_range = neigh.range();
    auto beg = neigh_range.begin;
    auto end = neigh_range.end;

    for (int_t d=0; d<num_dims_; ++d) {
      if (sector[d] != 0) {
        if (sector[d]<0) beg[d] -= num_ghost-1;
        end[d] = beg[d] + num_ghost;
      }
    }

    // rank and block id
    auto nblk = neigh.id();
    auto nrnk = owner(block_offsets, nblk);

    // added shared info
    for (int_t k=beg[2]; k!=end[2]; k++) {
      for (int_t j=beg[1]; j!=end[1]; j++) {
        for (int_t i=beg[0]; i!=end[0]; i++) {
          auto lid = cells.id(i, j, k);
          shared_info.emplace_back(lid, nrnk, nblk, lid);
        }
      }
    }

  } // sector neighbors
  
  //----------------------------------------------------------------------------
  // Finish intrablock connectiivty

  link_intrablock_cells();
  
  
}

////////////////////////////////////////////////////////////////////////////////
/// Build intrablock cell data structures
////////////////////////////////////////////////////////////////////////////////
void mesh_struct_block_t::link_intrablock_cells()
{
  intrablock_cell_to_cell_.clear();
  intrablock_cell_to_intrablock_faces_.clear();

  // make sure info exists
  auto it = cells2ghost_.find(num_dims_-1);
  if (it == cells2ghost_.end()) return;

  const auto & cells2ghost_face = it->second;

  auto & offsets = intrablock_cell_to_intrablock_faces_.offsets;
  auto & indices = intrablock_cell_to_intrablock_faces_.indices;
  
  auto nghost = cells2ghost_face.size();
  intrablock_cell_to_cell_.reserve(nghost);
  offsets.reserve(nghost+1);
  indices.reserve(nghost);

  offsets.emplace_back(0);

  for (auto & cm : cells2ghost_face) {
    auto cid = cm.first;
    const auto & fs = cm.second;
    intrablock_cell_to_cell_.emplace_back(cid);
    for (auto f : fs) indices.emplace_back(f);
    offsets.emplace_back(indices.size());
  }

  
}

////////////////////////////////////////////////////////////////////////////////
/// Vertices
////////////////////////////////////////////////////////////////////////////////
std::vector<real_t> mesh_struct_block_t::coordinates() const
{
  auto nvert = std::pow(2, num_dims_);
  std::vector<real_t> coords( nvert );

  auto nx = iterator_.ni;
  auto ny = iterator_.nj;
  auto nz = iterator_.nk;

  auto vertices = iterator_.vertices();

  if (num_dims_ == 1) {
    coords[0] = vertex_coords_[ vertices.id( 0, 0, 0) ];
    coords[1] = vertex_coords_[ vertices.id(nx, 0, 0) ];
  }
  else if (num_dims_ == 2) {
    int_t ids[] = {
      vertices.id( 0, 0, 0),
      vertices.id(nx, 0, 0),
      vertices.id(nx,ny, 0),
      vertices.id( 0,ny, 0) };
    for (int_t  i=0; i<nvert; ++i) {
      coords[i*2  ] = vertex_coords_[ids[i]*2];
      coords[i*2+1] = vertex_coords_[ids[i]*2 + 1];
    }
  }
  else if (num_dims_ == 3) {
    int_t ids[] = {
      vertices.id( 0, 0, 0),
      vertices.id(nx, 0, 0),
      vertices.id(nx,ny, 0),
      vertices.id( 0,ny, 0),
      vertices.id( 0, 0,nz),
      vertices.id(nx, 0,nz),
      vertices.id(nx,ny,nz),
      vertices.id( 0,ny,nz)
    };
    for (int_t  i=0; i<nvert; ++i) {
      coords[i*3  ] = vertex_coords_[ids[i]*3];
      coords[i*3+1] = vertex_coords_[ids[i]*3 + 1];
    }
  }

  return coords;
}

////////////////////////////////////////////////////////////////////////////////
/// Main block refinement
////////////////////////////////////////////////////////////////////////////////
std::vector<std::unique_ptr<mesh_struct_block_t>>
mesh_struct_block_t::refine(
    const amr_flags_t & flags,
    const range_ijk_t * ranges,
    int_t nrange) const
{
  std::vector<std::unique_ptr<mesh_struct_block_t>> new_blocks;
      
  int_t mult[] = {
    flags.has_refine(0) ? 2 : 1,
    flags.has_refine(1) ? 2 : 1,
    flags.has_refine(2) ? 2 : 1
  };

  auto vertices = iterator_.vertices();

  //----------------------------------------------------------------------------
  // One dimensional
  if (num_dims_ == 1) {

    assert(nrange == 2);
    UNUSED(nrange);

    auto vit = vertices.begin();
    auto v0 = vit.id();

    auto n0 = ranges[0].size(0)*2;
    auto n1 = ranges[1].size(0)*2;
    new_blocks.emplace_back(std::make_unique<mesh_struct_block_t>(region_, &n0, num_dims_));
    new_blocks.emplace_back(std::make_unique<mesh_struct_block_t>(region_, &n1, num_dims_));

    for (auto & b : new_blocks) {
      
      auto new_vit = b->iterator_.vertices().begin();
      b->vertex_coords_[new_vit.id()] = vertex_coords_[v0];

      for (int_t i=0; i<b->dim(0)/2; ++i) {
        vit.next();
        auto v1 = vit.id();
        auto xv0 = vertex_coords_[v0];
        auto xv1 = vertex_coords_[v1];
        auto xmid = 0.5*(xv0 + xv1);
        new_vit.next();
        b->vertex_coords_[new_vit.id()] = xmid;
        new_vit.next();
        b->vertex_coords_[new_vit.id()] = xv1;
        v0 = v1;
      }
    } // blocks

  }
  //----------------------------------------------------------------------------
  // Two dimensional
  else if (num_dims_ == 2) {

    assert(nrange == 2 || nrange == 4);
    UNUSED(nrange);

    for (int_t r=0; r<nrange; ++r) {
      const auto & rng = ranges[r];
      
      // create the new block
      long_t dims[] = { 
        rng.size(0) * mult[0], 
        rng.size(1) * mult[1],
      };
      auto b = std::make_unique<mesh_struct_block_t>(region_, dims, num_dims_);
      auto new_vertices = b->iterator_.vertices();

      for (int_t j=rng.begin[1], jj=0; j<rng.end[1]; ++j, jj+=mult[1]) {
        for (int_t i=rng.begin[0], ii=0; i<rng.end[0]; ++i, ii+=mult[0]) {

          auto v0 = vertices.id(i  , j  , 0);
          auto v1 = vertices.id(i+1, j  , 0);
          auto v2 = vertices.id(i+1, j+1, 0);
          auto v3 = vertices.id(i  , j+1, 0);

          real_t coefs[4][2];
          for (int_t d=0; d<2; ++d) {
            coefs[0][d] = vertex_coords_[ num_dims_*v0 + d ];
            coefs[1][d] = vertex_coords_[ num_dims_*v1 + d ] - coefs[0][d];
            coefs[2][d] = vertex_coords_[ num_dims_*v3 + d ] - coefs[0][d];
            coefs[3][d] = vertex_coords_[ num_dims_*v2 + d ] - coefs[0][d] - coefs[1][d] - coefs[2][d];
          }

          for (int_t jjj=0; jjj<=mult[1]; ++jjj) {
            auto yfact = static_cast<real_t>(jjj) / mult[1];
            for (int_t iii=0; iii<=mult[0]; ++iii) {
              auto xfact = static_cast<real_t>(iii) / mult[0];

              auto new_v = new_vertices.id(ii+iii, jj+jjj, 0);
              for (int_t d=0; d<2; ++d) {
                b->vertex_coords_[new_v * num_dims_ + d] =
                  coefs[0][d] +
                  coefs[1][d]*xfact +
                  coefs[2][d]*yfact +
                  coefs[3][d]*xfact*yfact;
              }

            } // ii
          } // jj


        } // i
      } // j

      new_blocks.emplace_back(std::move(b));
    } // blocks

  }
  //----------------------------------------------------------------------------
  // Three dimensional
  else if (num_dims_ == 3) {

    assert(nrange == 2 || nrange == 4 || nrange == 8);
    UNUSED(nrange);

    for (int_t r=0; r<nrange; ++r) {
      const auto & rng = ranges[r];
      
      // create the new block
      long_t dims[] = { 
        rng.size(0) * mult[0], 
        rng.size(1) * mult[1],
        rng.size(2) * mult[2]
      };
      auto b = std::make_unique<mesh_struct_block_t>(region_, dims, num_dims_);
      auto new_vertices = b->iterator_.vertices();

      for (int_t k=rng.begin[2], kk=0; k<rng.end[2]; ++k, kk+=mult[2]) {
        for (int_t j=rng.begin[1], jj=0; j<rng.end[1]; ++j, jj+=mult[1]) {
          for (int_t i=rng.begin[0], ii=0; i<rng.end[0]; ++i, ii+=mult[0]) {

            auto v0 = vertices.id(i  , j  , k  );
            auto v1 = vertices.id(i+1, j  , k  );
            auto v2 = vertices.id(i+1, j+1, k  );
            auto v3 = vertices.id(i  , j+1, k  );
            auto v4 = vertices.id(i  , j  , k+1);
            auto v5 = vertices.id(i+1, j  , k+1);
            auto v6 = vertices.id(i+1, j+1, k+1);
            auto v7 = vertices.id(i  , j+1, k+1);

            real_t coefs[8][3];
            for (int_t d=0; d<3; ++d) {
              coefs[0][d] = vertex_coords_[ num_dims_*v0 + d ];
              coefs[1][d] = vertex_coords_[ num_dims_*v1 + d ] - coefs[0][d];
              coefs[2][d] = vertex_coords_[ num_dims_*v3 + d ] - coefs[0][d];
              coefs[3][d] = vertex_coords_[ num_dims_*v4 + d ] - coefs[0][d];
              coefs[4][d] = vertex_coords_[ num_dims_*v2 + d ] - coefs[0][d] - coefs[1][d] - coefs[2][d];
              coefs[5][d] = vertex_coords_[ num_dims_*v5 + d ] - coefs[0][d] - coefs[1][d] - coefs[3][d];
              coefs[6][d] = vertex_coords_[ num_dims_*v7 + d ] - coefs[0][d] - coefs[2][d] - coefs[3][d];
              coefs[7][d] = vertex_coords_[ num_dims_*v6 + d ] - coefs[0][d] - coefs[1][d] - coefs[2][d]
                - coefs[3][d] - coefs[4][d] - coefs[5][d] - coefs[6][d];
            }

            for (int_t kkk=0; kkk<=mult[2]; ++kkk) {
              auto zfact = static_cast<real_t>(kkk) / mult[2];
              for (int_t jjj=0; jjj<=mult[1]; ++jjj) {
                auto yfact = static_cast<real_t>(jjj) / mult[1];
                for (int_t iii=0; iii<=mult[0]; ++iii) {
                  auto xfact = static_cast<real_t>(iii) / mult[0];

                  auto new_v = new_vertices.id(ii+iii, jj+jjj, kk+kkk);
                  for (int_t d=0; d<3; ++d) {
                    b->vertex_coords_[new_v * num_dims_ + d] =
                      coefs[0][d] +
                      coefs[1][d]*xfact +
                      coefs[2][d]*yfact +
                      coefs[3][d]*zfact +
                      coefs[4][d]*xfact*yfact +
                      coefs[5][d]*xfact*zfact + 
                      coefs[6][d]*yfact*zfact +
                      coefs[7][d]*xfact*yfact*zfact;
                  }

                } // ii
              } // jj
            } // kk


          } // i
        } // j
      } // k

      new_blocks.emplace_back(std::move(b));
    }

  } // dims    

  //----------------------------------------------------------------------------
  // Copy metadata

  for (int_t i=0; i<nrange; ++i) {
    const auto & rng = ranges[i];
    const auto & b = new_blocks[i];

    //--- Boundary conditions
    for (int_t d=0; d<num_dims_; ++d) {
      if (rng.begin[d] == 0) 
        b->add_boundary_mappings(0,d, boundary_ids(0,d));
      if (rng.end[d] == dims_[d]) 
        b->add_boundary_mappings(1,d, boundary_ids(1,d));
    }

    //--- copy metadata
    b->copy_metadata(*this);

  }

  return new_blocks;
}
  
////////////////////////////////////////////////////////////////////////////////
/// Main block coarsening
////////////////////////////////////////////////////////////////////////////////
std::unique_ptr<mesh_struct_block_t> mesh_struct_block_t::coarsen(
    int_t split,
    mesh_struct_block_t const * const * siblings)
{
  //----------------------------------------------------------------------------
  // Coarsen

  auto c0 = siblings[0];
  auto c1 = siblings[1];
  auto ndims = c0->num_dims();

  // compute dims
  std::vector<long_t> dims(ndims);
  for (int_t d=0; d<ndims; ++d)
    dims[d] = (d==split) ? (c0->dim(d) + c1->dim(d)) / 2 : c0->dim(d);
      
  auto new_block = std::make_unique<mesh_struct_block_t>(c0->region_, &dims[0], ndims);
  auto new_vertices = new_block->iterator_.vertices();

  int_t mult[3] = { 1, 1, 1 };
  mult[split] = 2;
 
  int_t start[] = {0, 0, 0};
  int_t end[] = {1, 1, 1};

  int_t nverts = std::pow(2, ndims);
  int_t new_v[8];
  int_t old_v[8];
   
  for (int_t c=0; c<2; ++c) {

    auto child = siblings[c];
    auto vertices = child->iterator_.vertices();

    for (int_t d=0; d<ndims; ++d)
      end[d] = start[d] + child->dim(d) / mult[d];

    for (int_t k=start[2], kk=0; k<end[2]; ++k, kk+=mult[2]) {
      for (int_t j=start[1], jj=0; j<end[1]; ++j, jj+=mult[1]) {
        for (int_t i=start[0], ii=0; i<end[0]; ++i, ii+=mult[0]) {

          new_v[0] = new_vertices.id(i  , j  , k  );
          new_v[1] = new_vertices.id(i+1, j  , k  );
          
          old_v[0] = vertices.id(ii        , jj  , kk  );
          old_v[1] = vertices.id(ii+mult[0], jj  , kk );

          for (int_t v=0; v<nverts; ++v)
            for (int_t d=0; d<ndims; ++d)
              new_block->vertex_coords_[new_v[v] * ndims + d] =
                child->vertex_coords_[old_v[v] * ndims + d];
        }
      }
    } // k
    
    start[split] = end[split];

  } // children
  
  //----------------------------------------------------------------------------
  // Copy metadata

  //--- Boundary conditions
  for (int_t d=0; d<ndims; ++d) {
    new_block->add_boundary_mappings(0,d, c0->boundary_ids(0,d));
    new_block->add_boundary_mappings(1,d, c0->boundary_ids(1,d));
  }

  //--- copy metadata
  new_block->copy_metadata(*c0);

  return new_block;

}

} // namespace
