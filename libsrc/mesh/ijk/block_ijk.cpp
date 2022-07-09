#include "block_ijk.hpp"
#include "mesh_ijk.hpp"

#include "utils/errors.hpp"

#include <algorithm>
#include <fstream>
#include <memory>

namespace prl {
  
using arr3_t = block_ijk_t::int3_t;

sector_map_t<arr3_t> block_ijk_t::sector_ctm_[3] = {
  {
    std::make_pair(sector_t{-1, 0, 0}, arr3_t{ 1, 0, 0}),
    std::make_pair(sector_t{ 1, 0, 0}, arr3_t{-1, 0, 0}),
  },
  {
    //--- faces
    std::make_pair(sector_t{-1, 0, 0}, arr3_t{ 1, 2, 0}),
    std::make_pair(sector_t{ 1, 0, 0}, arr3_t{-1, 2, 0}),
    std::make_pair(sector_t{ 0,-1, 0}, arr3_t{ 2, 1, 0}),
    std::make_pair(sector_t{ 0, 1, 0}, arr3_t{-2, 1, 0}),
    //--- corners 
    std::make_pair(sector_t{-1,-1, 0}, arr3_t{ 1, 2, 0}),
    std::make_pair(sector_t{ 1,-1, 0}, arr3_t{-1, 2, 0}),
    std::make_pair(sector_t{-1, 1, 0}, arr3_t{ 1,-2, 0}),
    std::make_pair(sector_t{ 1, 1, 0}, arr3_t{-1,-2, 0}),
  },
  {
    //--- faces
    std::make_pair(sector_t{-1, 0, 0}, arr3_t{ 1, 2, 3}),
    std::make_pair(sector_t{ 1, 0, 0}, arr3_t{-1, 3, 2}),
    std::make_pair(sector_t{ 0,-1, 0}, arr3_t{ 2, 3, 1}),
    std::make_pair(sector_t{ 0, 1, 0}, arr3_t{-2, 1, 3}),
    std::make_pair(sector_t{ 0, 0,-1}, arr3_t{ 3, 1, 2}),
    std::make_pair(sector_t{ 0, 0, 1}, arr3_t{-3, 2, 1}),
    //--- edges
    std::make_pair(sector_t{ 0,-1,-1}, arr3_t{ 1, 2, 3}),
    std::make_pair(sector_t{ 0, 1,-1}, arr3_t{ 1, 3,-2}),
    std::make_pair(sector_t{ 0,-1, 1}, arr3_t{ 1,-3, 2}),
    std::make_pair(sector_t{ 0, 1, 1}, arr3_t{ 1,-2,-3}),
    std::make_pair(sector_t{-1, 0,-1}, arr3_t{ 2, 3, 1}),
    std::make_pair(sector_t{ 1, 0,-1}, arr3_t{ 2,-1, 3}),
    std::make_pair(sector_t{-1, 0, 1}, arr3_t{ 2, 1,-3}),
    std::make_pair(sector_t{ 1, 0, 1}, arr3_t{ 2,-3,-1}),
    std::make_pair(sector_t{-1,-1, 0}, arr3_t{ 3, 1, 2}),
    std::make_pair(sector_t{ 1,-1, 0}, arr3_t{ 3, 2,-1}),
    std::make_pair(sector_t{-1, 1, 0}, arr3_t{ 3,-2, 1}),
    std::make_pair(sector_t{ 1, 1, 0}, arr3_t{ 3,-1,-2}),
    //--- corners 
    std::make_pair(sector_t{-1,-1,-1}, arr3_t{ 1, 2, 3}),
    std::make_pair(sector_t{ 1,-1,-1}, arr3_t{-1, 3, 2}),
    std::make_pair(sector_t{-1, 1,-1}, arr3_t{ 1, 3,-2}),
    std::make_pair(sector_t{ 1, 1,-1}, arr3_t{-1,-2, 3}),
    std::make_pair(sector_t{-1,-1, 1}, arr3_t{ 1,-3, 2}),
    std::make_pair(sector_t{ 1,-1, 1}, arr3_t{-1, 2,-3}),
    std::make_pair(sector_t{-1, 1, 1}, arr3_t{ 1,-2,-3}),
    std::make_pair(sector_t{ 1, 1, 1}, arr3_t{-1,-3,-2})

  }
};

////////////////////////////////////////////////////////////////////////////////
/// Determine the offset for a specific rotation
////////////////////////////////////////////////////////////////////////////////
int_t rotation2offset(int_t rot)
{ return rot < 0 ? -1 : 0; }

////////////////////////////////////////////////////////////////////////////////
/// Compute delta between two vertices
////////////////////////////////////////////////////////////////////////////////
void delta(int_t a, int_t b, const real_t * x, int_t ndims, real_t * dx)
{
  for (int_t d=0; d<ndims; ++d)
    dx[d] = x[a*ndims + d] - x[b*ndims + d];
}

////////////////////////////////////////////////////////////////////////////////
/// Compute delta between two vertices
////////////////////////////////////////////////////////////////////////////////
real_t set_dot(
    const real_t ** a,
    const real_t ** b,
    int_t irot,
    int_t ndims)
{
  std::vector<int_t> s2v(ndims);
  for (int_t i=0; i<ndims; ++i)
    s2v[i] = (i+ndims - irot) % ndims;
  
  real_t dot = 0;
  for (int_t i=0; i<ndims; ++i) {
    auto s2vi = s2v[i];
    for (int_t j=0; j<ndims; ++j) {
      dot += a[i][j]*b[s2vi][j];
    }
  }

  return dot;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for block class
////////////////////////////////////////////////////////////////////////////////
block_ijk_t::block_ijk_t(
    int_t id,
    const long_t * end, 
    const int_t * verts,
    size_t ndims,
    size_t nverts,
    mesh_ijk_t & mesh,
    const long_t * begin) :
  id_(id),
  label_(id),
  num_dims_(ndims),
  vertex_ids_(verts, verts+nverts),
  mesh_(mesh)
{

  std::copy_n(end, num_dims_, &range_.end[0]);
  if (begin) std::copy_n(begin, num_dims_, &range_.begin[0]);
  
  for (int_t d=0; d<num_dims_; ++d)
    dims_[d] = range_.end[d] - range_.begin[d];
    
  const auto coords = mesh_.coordinates();

  //--------------------------------------------------------------------------
  if (num_dims_ == 1) {
    
    // corners
    connect( -1, 0, 0, {verts[0]} );
    connect(  1, 0, 0, {verts[1]} );

    // direction vector
    for (int_t idir=-1; idir<2; idir+=2) {
      sector_t asec(idir, 0, 0), bsec(-idir, 0, 0);
      auto a = sector2verts_.at(asec).front();
      auto b = sector2verts_.at(bsec).front();
      block_ijk_vertex_t va(a, num_dims_);
      delta(b, a, &coords[0], num_dims_, va.direction(0));
      add_vertex(asec, std::move(va));
    }

  } // dims == 1
  //--------------------------------------------------------------------------
  else if (num_dims_ == 2) {
    
    // faces 
    connect(  0, -1, 0, {verts[0], verts[1]} );
    connect(  1,  0, 0, {verts[1], verts[2]} );
    connect(  0,  1, 0, {verts[3], verts[2]} );
    connect( -1,  0, 0, {verts[0], verts[3]} );

    // corners
    connect( -1, -1, 0, {verts[0]} );
    connect(  1, -1, 0, {verts[1]} );
    connect(  1,  1, 0, {verts[2]} );
    connect( -1,  1, 0, {verts[3]} );

    // direction vector
    for (int_t idir=-1; idir<2; idir+=2) {
      for (int_t jdir=-1; jdir<2; jdir+=2) {
        sector_t asec(idir, jdir, 0);
        auto a = sector2verts_.at(asec).front();
        block_ijk_vertex_t va(a, num_dims_);
        for (int_t d=0; d<2; ++d)
        {
          auto bsec = asec;
          bsec[d] *= -1;
          auto b = sector2verts_.at(bsec).front();
          delta(b, a, &coords[0], num_dims_, va.direction(d));
        }
        add_vertex(asec, std::move(va));
      }
    }

  } // dim==2
  //--------------------------------------------------------------------------
  else if (num_dims_ == 3) {

    // faces
    connect( -1,  0,  0, {verts[0], verts[4], verts[7], verts[3]} );
    connect(  1,  0,  0, {verts[1], verts[2], verts[6], verts[5]} );
    connect(  0, -1,  0, {verts[0], verts[1], verts[5], verts[4]} );
    connect(  0,  1,  0, {verts[3], verts[7], verts[6], verts[2]} );
    connect(  0,  0, -1, {verts[0], verts[3], verts[2], verts[1]} );
    connect(  0,  0,  1, {verts[4], verts[5], verts[6], verts[7]} );

    // edgees
    connect(  0, -1, -1, {verts[0], verts[1]} );
    connect(  0,  1, -1, {verts[3], verts[2]} );
    connect(  0, -1,  1, {verts[4], verts[5]} );
    connect(  0,  1,  1, {verts[7], verts[6]} );
    
    connect( -1,  0, -1, {verts[0], verts[3]} );
    connect(  1,  0, -1, {verts[1], verts[2]} );
    connect( -1,  0,  1, {verts[4], verts[7]} );
    connect(  1,  0,  1, {verts[5], verts[6]} );
    
    connect( -1, -1,  0, {verts[0], verts[4]} );
    connect(  1, -1,  0, {verts[1], verts[5]} );
    connect( -1,  1,  0, {verts[3], verts[7]} );
    connect(  1,  1,  0, {verts[2], verts[6]} );

    // corners
    connect( -1, -1, -1, {verts[0]} );
    connect(  1, -1, -1, {verts[1]} );
    connect(  1,  1, -1, {verts[2]} );
    connect( -1,  1, -1, {verts[3]} );
    connect( -1, -1,  1, {verts[4]} );
    connect(  1, -1,  1, {verts[5]} );
    connect(  1,  1,  1, {verts[6]} );
    connect( -1,  1,  1, {verts[7]} );

    // direction  vectors
    for (int_t idir=-1; idir<2; idir+=2) {
      for (int_t jdir=-1; jdir<2; jdir+=2) {
        for (int_t kdir=-1; kdir<2; kdir+=2) {
          sector_t asec(idir, jdir, kdir);
          auto a = sector2verts_.at(asec).front();
          block_ijk_vertex_t va(a, num_dims_);
          for (int_t d=0; d<3; ++d)
          {
            auto bsec = asec;
            bsec[d] *= -1;
            auto b = sector2verts_.at(bsec).front();
            delta(b, a, &coords[0], num_dims_, va.direction(d));
          }
          add_vertex(asec, std::move(va));
        }
      }
    }

  } // dims
  //--------------------------------------------------------------------------

}

////////////////////////////////////////////////////////////////////////////////
/// Constructor for block class
////////////////////////////////////////////////////////////////////////////////
block_ijk_t::block_ijk_t(
    block_ijk_t * parent,
    const range_ijk_t & range, 
    const int_t * verts,
    size_t ndims,
    size_t nverts,
    mesh_ijk_t & mesh) :
  block_ijk_t(-1, &range.end[0], verts, ndims, nverts, mesh, &range.begin[0])
{
  parent_ = parent;
  label_ = parent->label();
  
  level_ = parent->levels();
  level_[parent->split_dim()]++;

  spacing_ = parent->spacing_;
}

////////////////////////////////////////////////////////////////////////////////
/// Determine block range
////////////////////////////////////////////////////////////////////////////////
range_ijk_t block_ijk_t::range(const sector_t & sector) const
{
  range_ijk_t r;

  auto ndims = num_dims();
  for (int_t d=0; d<ndims; ++d) {
    if (sector[d] == 0) {
      r.begin[d] = 0;
      r.end[d] = dims_[d];
    }
    else if (sector[d] < 0) {
      r.begin[d] = 0;
      r.end[d] = 1;
    }
    else {
      r.begin[d] = dims_[d]-1;
      r.end[d] = dims_[d];
    }
  }

  return r;
}

range_ijk_t block_ijk_t::range_wrt_root(const sector_t & sector) const
{
  range_ijk_t r;

  auto ndims = num_dims();
  for (int_t d=0; d<ndims; ++d) {
    if (sector[d] == 0) {
      r.begin[d] = range_.begin[d];
      r.end[d] = range_.end[d];
    }
    else if (sector[d] < 0) {
      r.begin[d] = range_.begin[d];
      r.end[d] = range_.begin[d]+1;
    }
    else {
      r.begin[d] = range_.end[d]-1;
      r.end[d] = range_.end[d];
    }
  }

  return r;
}


////////////////////////////////////////////////////////////////////////////////
/// get coords
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::coords(real_t * xs) const
{
  const auto coord = mesh_.coordinates();
  auto ndims = num_dims();
  auto nverts = num_vertices();
  for (int_t v=0; v<nverts; ++v) {
    const auto xv = &coord[vertex_ids_[v]*ndims];
    for (int_t d=0; d<ndims; ++d)
        xs[v*ndims+d] = xv[d];
  }
}

std::vector<real_t> block_ijk_t::coords() const
{
  auto ndims = num_dims();
  auto nverts = num_vertices();
  std::vector<real_t> xv(ndims*nverts);
  coords(&xv[0]);
  return xv;
}
  
////////////////////////////////////////////////////////////////////////////////
/// Connnect a set of vertices to a sector
////////////////////////////////////////////////////////////////////////////////
result_t<sector_t> block_ijk_t::sector(std::vector<int_t> vs) const
{
  std::sort(vs.begin(), vs.end());
  auto it = verts2sector_.find(vs);
  if (it != verts2sector_.end())
    return {it->second, true};
  else
    return {false};
}
  
////////////////////////////////////////////////////////////////////////////////
/// Connnect a set of vertices to a sector
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::connect(int_t i, int_t j, int_t k, std::vector<int_t> conn)
{
  auto sect = sector_t{i, j, k};
  sector2verts_[sect] = conn;
  std::sort(conn.begin(), conn.end());
  verts2sector_[conn] = sect;
}
          
////////////////////////////////////////////////////////////////////////////////
/// Add a vertex
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::add_vertex(sector_t sec, block_ijk_vertex_t v) {
  vertices_.emplace(sec, std::move(v));
  vertex_map_.emplace(v.id_, sec);
}

////////////////////////////////////////////////////////////////////////////////
/// check if a certain neighbor exists
////////////////////////////////////////////////////////////////////////////////
bool block_ijk_t::has_neighbor_at(sector_t sec, const block_ijk_t & n) const
{ 
  auto it = neighbors_.find(sec);
  if (it != neighbors_.end()) {
    const auto & neighs = it->second;
    auto res = std::find_if(
        neighs.begin(),
        neighs.end(),
        [&](const auto & a) { return a.block()==n;} );
    return res != neighs.end();
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////////
/// Orient a block
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::orient(
    const std::vector<int_t> & sorted_vs,
    block_ijk_t & donor)
{

  // check if already added
  auto my_sector = verts2sector_.at(sorted_vs);

  if (has_neighbor_at(my_sector, donor)) return;

  auto sector_dims = my_sector.nnz();
  const auto & my_vs = sector2verts_.at(my_sector);
  auto vcommon = my_vs.front();
      
  auto ndims = num_dims();
  const auto & ref_map = sector_ctm_[ndims-1];

  auto donor_sector = donor.verts2sector_.at(sorted_vs);
  const auto & donor_vs = donor.sector2verts_.at(donor_sector);
  const auto & donor_dims = donor.dims_;
    
  // determine number of rotations
  auto it = std::find(donor_vs.begin(), donor_vs.end(), vcommon);
  if (it == donor_vs.end())
    THROW_ERROR("Could not find donor vertex despite supposed match.");
  auto irot = std::distance(donor_vs.begin(), it);

  // deteermine transformation between local and reference
  if (!ref_map.count(my_sector)) return;

  const auto & ctm_ref2me = ref_map.at(my_sector);
  transform_t ref2me( &ctm_ref2me[0], ndims );
  auto me2ref = transpose(ref2me);
  
  // flip the sector
  auto flip_sector = donor_sector;
  flip_sector.flip();
  
  // determine transformation between donor and reference
  auto ctm_ref2donor = ref_map.at(flip_sector);

  // outputs
  std::array<int_t, 3> my_offset = {0, 0, 0};
  std::array<int_t, 3> donor_offset = {0, 0, 0};
  std::unique_ptr<transform_t> me2donor;
  
  //----------------------------------------------------------------------------
  // Match 1d face
  //----------------------------------------------------------------------------
  if (ndims == 1) {
    
    // Adjustments based on rotation of the donor block
    constexpr int_t rotation_ctm[2] = { 1, -1};
 
    // determine transformation between donor and reference
    const auto & rotation = rotation_ctm[irot];
    ctm_rotate(&ctm_ref2donor[0], &rotation, ndims);
    transform_t ref2donor( &ctm_ref2donor[0], ndims );
    
    // relative transformation matrix
    me2donor = std::make_unique<transform_t>(ref2donor, me2ref);

    // offset for my block
    my_offset[0] = -std::max(0, my_sector[0]);
   
    // transform to local frame
    ref2me.transform(&my_offset[0]);

    // offset for donor block
    donor_offset[0] = std::max(0, donor_sector[0]);
  
    // transform donor offset to donor frame
    ref2donor.transform(&donor_offset[0]);
  }
  //----------------------------------------------------------------------------
  // Match a corner
  //----------------------------------------------------------------------------
  else if (sector_dims == ndims) {
  
    // Adjustents based on rotation of the donor block
    constexpr int_t rotation_ctm[3][3][3] = {
      { { 0, 0, 0}, { 0, 0, 0}, { 0, 0, 0} },
      { { 1, 2, 0}, { 2, 1, 0}, { 0, 0, 0} },
      { { 1, 2, 3}, { 3, 1, 2}, { 2, 3, 1} }
    };

    // Vectors in the reference i, j, and k directions for the
    // local block.
    const auto & my_dirs = vertices_.at(my_sector).dirs_;
    const real_t * my_vecs[3];
    for (int_t i=0; i<ndims; ++i) 
      my_vecs[i] = &(my_dirs[std::abs(ctm_ref2me[i]) - 1])[0];

    // Vectors in the reference i, j, and k directions for the
    // donor block
    const auto & donor_dirs = donor.vertices_.at(donor_sector).dirs_;
    const real_t * donor_vecs[3];
    for (int_t i=0; i<ndims; ++i)
      donor_vecs[i] = &(donor_dirs[std::abs(ctm_ref2donor[i]) - 1])[0];

    // Weights of a dot-product between different orientations of
    // the vectors.
    real_t weights[3];
    for (int_t i=0; i<ndims; ++i)
      weights[i] = set_dot(my_vecs, donor_vecs, i, ndims);

    // The minimum weight gives the closest alignment
    auto it = std::min_element(weights, weights+ndims);
    auto irot = std::distance(weights, it);


    // determine transformation between donor and reference
    const auto & rotation = rotation_ctm[ndims-1][irot];
    ctm_rotate(&ctm_ref2donor[0], &rotation[0], ndims);
    transform_t ref2donor( &ctm_ref2donor[0], ndims );
    
    // relative transformation matrix
    me2donor = std::make_unique<transform_t>(ref2donor, me2ref);

    // offset for my block
    for (int_t i=0; i<ndims; ++i)
      my_offset[i] = std::max(my_sector[i], 0);

    // offset for donor block
    for (int_t i=0; i<ndims; ++i)
      donor_offset[i] = std::max(donor_sector[i], 0);
  }
  //----------------------------------------------------------------------------
  // Match an edge in 3d
  //----------------------------------------------------------------------------
  else if (sector_dims == 2) {

    // Adjustents based on rotation of the donor block
    constexpr int_t rotation_ctm[2][3] = { { 1, 2, 3}, {-1, 3, 2} };

    // determine transformation between donor and reference
    const auto & rotation = rotation_ctm[irot];
    ctm_rotate(&ctm_ref2donor[0], &rotation[0], ndims);
    transform_t ref2donor( &ctm_ref2donor[0], ndims );
    
    // relative transformation matrix
    me2donor = std::make_unique<transform_t>(ref2donor, me2ref);

    // offset for my block
    my_offset[0] = rotation2offset(ctm_ref2me[0]);
    my_offset[1] = rotation2offset(ctm_ref2me[1]);
    my_offset[2] = rotation2offset(ctm_ref2me[2]);
   
    // transform to local frame
    ref2me.transform(&my_offset[0]);

    // offset for donor block
    const auto & ctm_donor2ref = ref_map.at(donor_sector);
    if (irot == 0) {
      donor_offset[0] = 0; 
      donor_offset[1] = -rotation2offset(ctm_donor2ref[1]);
      donor_offset[2] = -rotation2offset(ctm_donor2ref[2]);
    }
    else {
      donor_offset[0] = -1;
      donor_offset[1] = -rotation2offset(ctm_donor2ref[2]);
      donor_offset[2] = -rotation2offset(ctm_donor2ref[1]);
    }
    
    // transform donor offset to donor frame
    ref2donor.transform(&donor_offset[0]);

  }
  //----------------------------------------------------------------------------
  // Match a face in 2d
  //----------------------------------------------------------------------------
  else if (ndims == 2 && sector_dims == 1) {
    
    // Adjustments based on rotation of the donor block
    constexpr int_t rotation_ctm[2][2] = { { 1, 2}, {1, -2} };
 
    // determine transformation between donor and reference
    const auto & rotation = rotation_ctm[irot];
    ctm_rotate(&ctm_ref2donor[0], &rotation[0], ndims);
    transform_t ref2donor( &ctm_ref2donor[0], ndims );
    
    // relative transformation matrix
    me2donor = std::make_unique<transform_t>(ref2donor, me2ref);

    // offset for my block
    my_offset[0] = -(std::max(0, my_sector[0]) + std::max(0, my_sector[1]));
    my_offset[1] = 0;
   
    // transform to local frame
    ref2me.transform(&my_offset[0]);
    
    // offset caused by rotation of fdonor block
    int_t rotation_offset[2] = {
      rotation2offset(rotation[0]),
      rotation2offset(rotation[1])
    };

    // offset for donor block
    donor_offset[0] = std::max(0, donor_sector[0]) + std::max(0, donor_sector[1]);
    donor_offset[1] = rotation_offset[1];
  
    // transform donor offset to donor frame
    ref2donor.transform(&donor_offset[0]);
  }
  //----------------------------------------------------------------------------
  // Match a face in 3d
  //----------------------------------------------------------------------------
  else if (ndims == 3 && sector_dims == 1) {
  
    // Adjustments based on rotation of the donor block
    constexpr int_t rotation_ctm[4][3] = { {1, 2, 3}, {1, 3, -2}, {1, -2, -3}, {1, -3, 2} };
 
    // determine transformation between donor and reference
    const auto & rotation = rotation_ctm[irot];
    ctm_rotate(&ctm_ref2donor[0], &rotation[0], ndims);
    transform_t ref2donor( &ctm_ref2donor[0], ndims );
    
    // relative transformation matrix
    me2donor = std::make_unique<transform_t>(ref2donor, me2ref);

    // offset for my block
    my_offset[0] = 
      -(std::max(0, my_sector[0]) + std::max(0, my_sector[1]) + std::max(0, my_sector[2]));
    my_offset[1] = 0;
    my_offset[2] = 0;
   

    // transform to local frame
    ref2me.transform(&my_offset[0]);
    
    // offset caused by rotation of fdonor block
    int_t rotation_offset[3] = {
      rotation2offset(rotation[0]),
      rotation2offset(rotation[1]),
      rotation2offset(rotation[2])
    };

    // offset for donor block
    donor_offset[0] =
      std::max(0, donor_sector[0]) +
      std::max(0, donor_sector[1]) +
      std::max(0, donor_sector[2]);
    donor_offset[1] = rotation_offset[1];
    donor_offset[2] = rotation_offset[2];
  
    // transform donor offset to donor frame
    ref2donor.transform(&donor_offset[0]);
  }
  
  //----------------------------------------------------------------------------
  // Finish up
  //----------------------------------------------------------------------------
  if (me2donor) {

    auto me2donor_verts = *me2donor;
    auto my_offset_verts = my_offset;
    auto donor_offset_verts = donor_offset;

    //--- Local offsets
    // multiply by dims
    for (int_t d=0; d<ndims; ++d) {
      my_offset      [d] *= dims_[d]-1;
      my_offset_verts[d] *= dims_[d];
    }
      
    // put in donor frame
    me2donor->transform(&my_offset[0]);
    me2donor->transform(&my_offset_verts[0]);

    //--- Donor offsets
    for (int_t d=0; d<ndims; ++d) {
      // mulitply by dims
      donor_offset      [d] *= donor_dims[d]-1;
      donor_offset_verts[d] *= donor_dims[d];
      // add jump accross boundary
      donor_offset[d] += donor_sector[d];
      // get relative ofset
      me2donor     ->offset(d) = donor_offset      [d] - my_offset      [d];
      me2donor_verts.offset(d) = donor_offset_verts[d] - my_offset_verts[d];
    }
    
    //--- Compute ranges
    auto my_range = range(my_sector);
    auto donor_range = donor.range(donor_sector); 

    // set neighbor info
    add_neighbor(
        my_sector,
        donor,
        *me2donor,
        me2donor_verts,
        my_range);
    me2donor->rotate(&my_sector[0]);
    donor.add_shared(*this, my_sector, donor_range);

    me2donor->reverse();
    me2donor_verts.reverse();

    donor.add_neighbor(
        donor_sector,
        *this,
        *me2donor,
        me2donor_verts,
        donor_range);
    me2donor->rotate(&donor_sector[0]);
    add_shared(donor, donor_sector, my_range);

  } // me2donor exists
}

////////////////////////////////////////////////////////////////////////////////
/// Root block constructor
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::set_linear_spacing()
{
  auto ndims = num_dims();

  //--- One dimensions
  if (ndims == 1) {
    spacing_ = [](const auto xv, const auto i, const auto n, auto x) 
    { 
      auto lower = xv[0];
      auto delta = xv[1] - lower;
      x[0] = lower + i[0] * delta / n[0];
    };
  }
  //--- Two dimensions
  else if (ndims == 2) {
    spacing_ = [](const auto xv, const auto i, const auto n, auto x) 
    { 
      auto xa = xv[ 2*0 + 0];  auto ya = xv[ 2*0 + 1];
      auto xb = xv[ 2*1 + 0];  auto yb = xv[ 2*1 + 1];
      auto xc = xv[ 2*2 + 0];  auto yc = xv[ 2*2 + 1];
      auto xd = xv[ 2*3 + 0];  auto yd = xv[ 2*3 + 1];
      auto xfact = static_cast<real_t>(i[0]) / n[0];
      auto yfact = static_cast<real_t>(i[1]) / n[1];
      auto x0 = (1-xfact)*xa + xfact*xb;
      auto x1 = (1-xfact)*xd + xfact*xc;
      auto y0 = (1-yfact)*ya + yfact*yd;
      auto y1 = (1-yfact)*yb + yfact*yc;
      x[0] = (1-yfact)*x0 + yfact*x1;
      x[1] = (1-xfact)*y0 + xfact*y1;
    };
  }
  //--- Three dimensions
  else if (ndims == 3) {
    spacing_ = [](const auto xv, const auto i, const auto n, auto x)
    {
      // trilinear mapping
      real_t coefs[8][3];
      for (int_t d=0; d<3; ++d) {
        coefs[0][d] = xv[ 3*0 + d ];
        coefs[1][d] = xv[ 3*1 + d ] - coefs[0][d];
        coefs[2][d] = xv[ 3*3 + d ] - coefs[0][d];
        coefs[3][d] = xv[ 3*4 + d ] - coefs[0][d];
        coefs[4][d] = xv[ 3*2 + d ] - coefs[0][d] - coefs[1][d] - coefs[2][d];
        coefs[5][d] = xv[ 3*5 + d ] - coefs[0][d] - coefs[1][d] - coefs[3][d];
        coefs[6][d] = xv[ 3*7 + d ] - coefs[0][d] - coefs[2][d] - coefs[3][d];
        coefs[7][d] = xv[ 3*6 + d ] - coefs[0][d] - coefs[1][d] - coefs[2][d]
          - coefs[3][d] - coefs[4][d] - coefs[5][d] - coefs[6][d];
      }
      auto xfact = static_cast<real_t>(i[0]) / n[0];
      auto yfact = static_cast<real_t>(i[1]) / n[1];
      auto zfact = static_cast<real_t>(i[2]) / n[2];
      for (int_t d=0; d<3; ++d) {
        x[d] =
          coefs[0][d] +
          coefs[1][d]*xfact +
          coefs[2][d]*yfact +
          coefs[3][d]*zfact +
          coefs[4][d]*xfact*yfact +
          coefs[5][d]*xfact*zfact + 
          coefs[6][d]*yfact*zfact +
          coefs[7][d]*xfact*yfact*zfact;
      }
    };
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Split a block
////////////////////////////////////////////////////////////////////////////////
bool block_ijk_t::split_block(size_t min) 
{

  //--- At leaf
  if (children_.empty()) {

    min = std::max<size_t>(min, 1);
    auto it = std::max_element(&dims_[0], &dims_[num_dims_]);
    split_ = std::distance(&dims_[0], it);
    
    if (static_cast<size_t>(dims_[split_]) < 2*min) return true;
    
    long_t splt = static_cast<real_t>(dims_[split_]) / 2 + 0.5;

    range_ijk_t ranges[2];

    {
      auto tmp = range_.end;
      tmp[split_] = range_.begin[split_] + splt;
      ranges[0] = range_ijk_t(num_dims_, &range_.begin[0], &tmp[0]);
    }
    {
      auto tmp = range_.begin;
      tmp[split_] += splt;
      ranges[1] = range_ijk_t(num_dims_, &tmp[0], &range_.end[0]);
    }

    split_block(ranges, 2);
    
    return false;
  }
  //--- Continue desecnnding until we reach a leaf
  else {

    bool ret = true;
    for (auto & c : children_)
      ret &= c->split_block(min);
    return ret;

  }

}


////////////////////////////////////////////////////////////////////////////////
/// void Split a block
////////////////////////////////////////////////////////////////////////////////
//std::vector<std::unique_ptr<block_ijk_t>>
void block_ijk_t::split_block(const range_ijk_t * ranges, int_t nranges)
{

  // use root spacing function (block ranges are relative to the root)
  auto rt = root();
  auto root_coords = rt->coords();
      
  children_.reserve(children_.size() + nranges);

  const auto & begin = range_.begin;
  const auto & end = range_.end;
  
  //----------------------------------------------------------------------------
  // Make a map to search for corners
  using long3_t = std::array<long_t,3>;

  struct hash_t {
    long3_t end;
    hash_t(const long3_t & e) : end(e) {}
    size_t operator()(const long3_t & point) const
    { return  point[0] + (end[0]+1)*(point[1] + (end[1]+1)*point[2] + point[1]); }
  };

  size_t bucket_count = nranges * std::pow(2, num_dims_); 
  std::unordered_map< long3_t, int_t, hash_t > vertex_map(bucket_count, hash_t{end}); 
  
  //----------------------------------------------------------------------------
  
  auto create_vertex = [&](const auto & pos)
  { 
    auto apos = long3_t{pos[0], pos[1], pos[2]};
    auto it = vertex_map.find(apos);
    if (it == vertex_map.end()) {
      real_t x[3];
      rt->spacing(&root_coords[0], &pos[0], rt->dims(), x);
      auto v = mesh_.add_vertex(x);
      vertex_map.emplace(apos, v);
      return v;
    }
    else {
      return it->second;
    }
  };
  
  auto create_block = [&](const auto & range, const auto & verts, const auto nverts)
  {
    auto b = std::make_unique<block_ijk_t>(this, range, verts, num_dims_, nverts, mesh_);
    children_.emplace_back( std::move(b) );
  };
  
  //----------------------------------------------------------------------------
  if (num_dims_ == 1) {
  
    vertex_map.emplace( long3_t{begin[0], 0, 0}, vertex_ids_[0]);
    vertex_map.emplace( long3_t{  end[0], 0, 0}, vertex_ids_[1]);
    
    int_t verts[2];
    
    for (int_t i=0; i<nranges; ++i) {
      const auto & r = ranges[i];
      verts[0] = create_vertex(&r.begin[0]);
      verts[1] = create_vertex(&r.end[0]);
      create_block(r, &verts[0], 2);
    }
    
  }
  //----------------------------------------------------------------------------
  else if (num_dims_ == 2) {

    vertex_map.emplace( long3_t{begin[0], begin[1], 0}, vertex_ids_[0]);
    vertex_map.emplace( long3_t{  end[0], begin[1], 0}, vertex_ids_[1]);
    vertex_map.emplace( long3_t{  end[0],   end[1], 0}, vertex_ids_[2]);
    vertex_map.emplace( long3_t{begin[0],   end[1], 0}, vertex_ids_[3]);
    
    int_t verts[4];

    for (int_t i=0; i<nranges; ++i) {
      const auto & r = ranges[i];
      auto tmp = r.begin;
      verts[0] = create_vertex(&tmp[0]);
      tmp[0] = r.end[0];
      verts[1] = create_vertex(&tmp[0]);
      tmp[1] = r.end[1];
      verts[2] = create_vertex(&tmp[0]);
      tmp[0] = r.begin[0];
      verts[3] = create_vertex(&tmp[0]);
      create_block(r, &verts[0], 4);
    }
  
  }
  //----------------------------------------------------------------------------
  else {

    vertex_map.emplace( long3_t{begin[0], begin[1], begin[2]}, vertex_ids_[0]);
    vertex_map.emplace( long3_t{  end[0], begin[1], begin[2]}, vertex_ids_[1]);
    vertex_map.emplace( long3_t{  end[0],   end[1], begin[2]}, vertex_ids_[2]);
    vertex_map.emplace( long3_t{begin[0],   end[1], begin[2]}, vertex_ids_[3]);
    vertex_map.emplace( long3_t{begin[0], begin[1],   end[2]}, vertex_ids_[4]);
    vertex_map.emplace( long3_t{  end[0], begin[1],   end[2]}, vertex_ids_[5]);
    vertex_map.emplace( long3_t{  end[0],   end[1],   end[2]}, vertex_ids_[6]);
    vertex_map.emplace( long3_t{begin[0],   end[1],   end[2]}, vertex_ids_[7]);
    
    int_t verts[8];

    for (int_t i=0; i<nranges; ++i) {
      const auto & r = ranges[i];
      // -ve z
      auto tmp = r.begin;
      verts[0] = create_vertex(&tmp[0]);
      tmp[0] = r.end[0];
      verts[1] = create_vertex(&tmp[0]);
      tmp[1] = r.end[1];
      verts[2] = create_vertex(&tmp[0]);
      tmp[0] = r.begin[0];
      verts[3] = create_vertex(&tmp[0]);
      // +ve z
      tmp = r.begin;
      tmp[2] = r.end[2];
      verts[4] = create_vertex(&tmp[0]);
      tmp[0] = r.end[0];
      verts[5] = create_vertex(&tmp[0]);
      tmp[1] = r.end[1];
      verts[6] = create_vertex(&tmp[0]);
      tmp[0] = r.begin[0];
      verts[7] = create_vertex(&tmp[0]);
      // create block
      create_block(r, &verts[0], 8);
    }
  
  } // dims
}

////////////////////////////////////////////////////////////////////////////////
/// Make a range relative to me
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::make_relative(range_ijk_t & rng) const
{
  for (int_t dim=0; dim<num_dims_; ++dim) {
    rng.begin[dim] -= range_.begin[dim];
    rng.end[dim]   -= range_.begin[dim];
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Add a block neighbor
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::add_internal_neighbor(
    block_ijk_t & neigh,
    const sector_t & my_sector,
    range_ijk_t my_range)
{
  
  //--- Compute a known position in target frame
  std::array<long_t,3> my_cell_pos = {0, 0, 0}, my_vert_pos = {0, 0, 0};
  for (int_t d=0; d<num_dims_; ++d) {
    my_cell_pos[d] = my_range.begin[d] + my_sector[d];
    my_vert_pos[d] = my_sector[d]>0 ? my_range.begin[d]+1 : my_range.begin[d];
  }

  //--- already in donor frame
  auto neigh_cell_pos = my_cell_pos;
  auto neigh_vert_pos = my_vert_pos;

  //--- make positions relative
  make_relative(&my_cell_pos[0]);
  make_relative(&my_vert_pos[0]);

  neigh.make_relative(&neigh_cell_pos[0]);
  neigh.make_relative(&neigh_vert_pos[0]);

  //--- compute transformation
  ctm_t me2neigh_cells(num_dims_), me2neigh_verts(num_dims_);
  
  //--- compute offsets
  for (int_t d=0; d<num_dims_; ++d) {
    me2neigh_cells.offset(d) = neigh_cell_pos[d] - my_cell_pos[d];
    me2neigh_verts.offset(d) = neigh_vert_pos[d] - my_vert_pos[d];
  }

  //--- Make ranges relative
  make_relative(my_range);
 
  //--- Add my neighbor 
  add_neighbor(my_sector, neigh, me2neigh_cells, me2neigh_verts, my_range);

  //--- Add shared info
  for (int_t d=0; d<num_dims_; ++d) {
    my_range.begin[d] += my_sector[d];
    my_range.end[d] += my_sector[d];
  }
  me2neigh_cells.transform(my_range);

  neigh.add_shared(*this, my_sector, my_range);
}

////////////////////////////////////////////////////////////////////////////////
/// Add a block neighbor
/// \note neigh_range is in the donor frame
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::add_external_neighbor(
    block_ijk_t & neigh,
    sector_t my_sector,
    range_ijk_t my_range,
    transform_t me2neigh_cells,
    transform_t me2neigh_verts)
{

  //--- Compute a known position in target frame
  std::array<long_t,3> my_cell_pos, my_vert_pos;
  for (int_t d=0; d<num_dims_; ++d) {
    my_cell_pos[d] = my_range.begin[d] + my_sector[d];
    my_vert_pos[d] = my_sector[d]>0 ? my_range.begin[d]+1 : my_range.begin[d];
  }

  //--- transform it to donor frame
  auto neigh_cell_pos = my_cell_pos;
  auto neigh_vert_pos = my_vert_pos;
  me2neigh_cells.transform(&neigh_cell_pos[0]);
  me2neigh_verts.transform(&neigh_vert_pos[0]);

  //--- make positions relative
  make_relative(&my_cell_pos[0]);
  make_relative(&my_vert_pos[0]);

  neigh.make_relative(&neigh_cell_pos[0]);
  neigh.make_relative(&neigh_vert_pos[0]);
  
  //--- Make initial transformation
  me2neigh_cells.zero_offsets();
  me2neigh_verts.zero_offsets();

  auto new_neigh_cell_pos = my_cell_pos;
  auto new_neigh_vert_pos = my_vert_pos;

  me2neigh_cells.transform(&new_neigh_cell_pos[0]);
  me2neigh_verts.transform(&new_neigh_vert_pos[0]);

  //--- compute offsets
  for (int_t d=0; d<num_dims_; ++d) {
    me2neigh_cells.offset(d) = neigh_cell_pos[d] - new_neigh_cell_pos[d];
    me2neigh_verts.offset(d) = neigh_vert_pos[d] - new_neigh_vert_pos[d];
  }

  //--- Make ranges relative
  make_relative(my_range);

  //--- Add targets neighbor
  add_neighbor(
      my_sector,
      neigh,
      me2neigh_cells,
      me2neigh_verts,
      my_range);
  
  //--- Add shared info
  for (int_t d=0; d<num_dims_; ++d) {
    my_range.begin[d] += my_sector[d];
    my_range.end[d] += my_sector[d];
  }
  me2neigh_cells.rotate(&my_sector[0]);
  me2neigh_cells.transform(my_range);
  neigh.add_shared(*this, my_sector, my_range);
}


////////////////////////////////////////////////////////////////////////////////
/// Number leaves
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::number(int_t & counter)
{
  if (children_.empty()) {
    id_ = counter++;
  }
  else {
    id_ = -1;
    for (const auto & c : children_)
      c->number(counter);
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Get all the leaves
////////////////////////////////////////////////////////////////////////////////
void block_ijk_t::get_leaves(std::vector<const block_ijk_t*> & list) const
{
  descend([&](const auto & node) {
    list.emplace_back(&node);
  });
}

std::vector<const block_ijk_t*> block_ijk_t::get_leaves() const
{
  auto n = count_leaves();
  std::vector<const block_ijk_t*> leaves;
  leaves.reserve(n);
  get_leaves(leaves);
  return leaves;
}

////////////////////////////////////////////////////////////////////////////////
/// Count leaves
////////////////////////////////////////////////////////////////////////////////
size_t block_ijk_t::count_leaves() const
{
  size_t n=0;
  descend( [&](const auto &){ n++; } );
  return n;
}

////////////////////////////////////////////////////////////////////////////////
/// Find the largest block
////////////////////////////////////////////////////////////////////////////////
const block_ijk_t * block_ijk_t::find_largest_block() const
{
  size_t size = 0;
  const block_ijk_t * largest = nullptr;

  descend( [&](const auto & node) {
    auto nsize = node.size();
    if (nsize > size) {
      largest = &node;
      size = nsize;
    }
  });

  return largest;

}

////////////////////////////////////////////////////////////////////////////////
/// Find children
////////////////////////////////////////////////////////////////////////////////
result_t<block_ijk_t::child_list_t::const_iterator>
block_ijk_t::find_child(const block_ijk_t * c) const
{
  auto it = std::find_if(
      children_.begin(),
      children_.end(),
      [=](const auto & child) { return child.get() == c; });
  if (it == children_.end())
    return {false};
  else
    return {it, true};
}

////////////////////////////////////////////////////////////////////////////////
/// Find left and right nodes
////////////////////////////////////////////////////////////////////////////////
const block_ijk_t * block_ijk_t::left(const block_ijk_t * c) const
{
  auto res = find_child(c);
  if (res && res.value != children_.begin())
    return std::prev(res.value)->get();
  else
    return nullptr;
}

const block_ijk_t * block_ijk_t::right(const block_ijk_t * c) const
{
  auto res = find_child(c);
  if (res && res.value != std::prev(children_.end()))
    return std::next(res.value)->get();
  else
    return nullptr;
}
  
////////////////////////////////////////////////////////////////////////////////
/// Find blocks belonging to a sector
////////////////////////////////////////////////////////////////////////////////
std::vector<const block_ijk_t *>
block_ijk_t::find_sectors(const sector_t & sector) const
{
  std::vector<const block_ijk_t *> list;
  find_sectors(sector, list);
  return list;
}

void block_ijk_t::find_sectors(
    const sector_t & sector,
    std::vector<const block_ijk_t *> & list) const
{ _find_sectors(this, sector, list); }

void block_ijk_t::_find_sectors(
    const block_ijk_t * root,
    const sector_t & sector,
    std::vector<const block_ijk_t *> & list) const
{
  if (children_.empty()) {
    list.emplace_back(this);
  }
  else {
    // search dim is same as split
    auto split = sector[split_];
    if (split) {
      const auto & child = split<0 ? children_.front() : children_.back();
      child->_find_sectors(root, sector, list);
    }
    // split along another dim
    else {
      for (const auto & c : children_) {
        c->_find_sectors(root, sector, list);
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
/// Find the neighbor blocks
////////////////////////////////////////////////////////////////////////////////
std::vector<const block_ijk_t *>
block_ijk_t::find_neighbors(const sector_t & sector) const
{
  std::vector<const block_ijk_t *> list;
  find_neighbors(sector, list);
  return list;
}

void block_ijk_t::find_neighbors(
    const sector_t & sector,
    std::vector<const block_ijk_t *> & list) const
{
  std::array<bool,3> mask;
  sector.copy(&mask[0]);
  _find_neighbors(this, sector, list, mask);
}


void block_ijk_t::_find_neighbors(
    const block_ijk_t * root,
    const sector_t & sector,
    std::vector<const block_ijk_t *> & list,
    std::array<bool,3> & mask) const
{

  auto has_left = std::accumulate(&mask[0], &mask[3], 0);

  // go up until the split matches 
  if (parent_ && has_left) {

    // if the split matches, mark it
    auto split_dim = parent_->split_;
    auto split_mask = sector[split_dim];
    if (split_mask && mask[split_dim]) {
        
      // if there is a neighbor to cross, cross it
      auto neighbor = (split_mask<0) ? parent_->left(this) : parent_->right(this);
      if (neighbor) {

        mask[split_dim] = false;
        has_left = std::accumulate(&mask[0], &mask[3], 0);

        if (!has_left) {
          neighbor->find_neighbors_descend(root->range_, sector, list);
          return;
        } // none left
        
      } // neighbor

    } // split

    // keep going up if we have more neighbors to find
    parent_->_find_neighbors(root, sector, list, mask);

  }
}

std::vector<const block_ijk_t *>
block_ijk_t::find_neighbors_descend(const sector_t & sector) const
{
  std::vector<const block_ijk_t *> list;
  find_neighbors_descend(sector, list);
  return list;
}

void block_ijk_t::find_neighbors_descend(
    const sector_t & sector,
    std::vector<const block_ijk_t *> & list) const
{ find_neighbors_descend(range_, sector, list); }

std::vector<const block_ijk_t *>
block_ijk_t::find_neighbors_descend(
    const range_ijk_t & range,
    const sector_t & sector) const
{
  std::vector<const block_ijk_t *> list;
  find_neighbors_descend(range, sector, list);
  return list;
}

void block_ijk_t::find_neighbors_descend(
    const range_ijk_t & range,
    const sector_t & sector,
    std::vector<const block_ijk_t *> & list) const
{
  if (children_.empty()) {
    if (do_ranges_intersect(num_dims_, sector, range, range_))
      list.emplace_back(this);
  }
  else {
    //--- search dim is same as split, find the closest block to it,
    //--- and stay to one side of the cut plane
    auto split = sector[split_];
    if (split) {
      if (split<0) {
        auto split_plane = range.begin[split_];
        for (const auto & c : children_) {
          if (c->begin(split_) < split_plane && c->end(split_) >= split_plane) {
            c->find_neighbors_descend(range, sector, list);
            return;
          }
        }
      } // -1
      else {
        auto split_plane = range.end[split_];
        for (const auto & c : children_) {
          if (c->begin(split_) <= split_plane && c->end(split_) > split_plane) {
            c->find_neighbors_descend(range, sector, list);
            return;
          }
        }
      } // 1
    }
    // split along another dim
    else {
      for (const auto & c : children_) {
        c->find_neighbors_descend(range, sector, list);
      }
    }
  }
}


  
} // namespace
