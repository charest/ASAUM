#ifndef LOGICAL_ITERATORS_HPP
#define LOGICAL_ITERATORS_HPP

#include "config.hpp"

namespace prl {

struct logical_vertex_t;
struct vertex_iterator_t;

struct logical_edge_t;
struct edge_iterator_t;

struct logical_face_t;
struct face_iterator_t;

struct logical_cell_t;
struct cell_iterator_t;

////////////////////////////////////////////////////////////////////////////////
/// \brief The main logical structured block iterator.
///
/// It is used to produce cell/face iterators
////////////////////////////////////////////////////////////////////////////////
struct logical_block_t {

  int_t nd = 0,
        ni = 0,
        nj = 0,
        nk = 0;

  logical_block_t() = default;

  logical_block_t(int_t d, int_t x, int_t y, int_t z) :
    nd(d), ni(x), nj(y), nk(z)
  {}
  
  logical_block_t(int_t d, const int_t * x) :
    nd(d), ni(x[0]), nj(d>1 ? x[1] : 1), nk(d>2 ? x[2] : 1)
  {}
  
  logical_block_t(int_t d, const long_t * x) :
    nd(d), ni(x[0]), nj(d>1 ? x[1] : 1), nk(d>2 ? x[2] : 1)
  {}
  logical_vertex_t vertices() const;
  
  logical_edge_t edges(int_t i) const;

  logical_face_t faces(int_t d) const;

  logical_cell_t cells() const;
};

////////////////////////////////////////////////////////////////////////////////
/// Class to iterate over vertices
////////////////////////////////////////////////////////////////////////////////
struct logical_vertex_t {
  const logical_block_t & block;
  
  int_t ni=0, nj=0, nk=0;
  int_t strides[3] = {0, 0, 0};

  logical_vertex_t(const logical_block_t & blk) :
    block(blk), strides{1, block.ni+1, (block.ni+1)*(block.nj+1)}
  {
    auto nd = blk.nd;
    ni = blk.ni+1;
    if (nd>1) nj = blk.nj+1;
    if (nd>2) nk = blk.nk+1;
  }
  
  int_t id(int_t i, int_t j, int_t k) const
  { return i*strides[0] + j*strides[1] + k*strides[2]; }
  
  vertex_iterator_t begin() const;

  int_t size() const
  { 
    if (block.nd == 1)
      return block.ni+1;
    else if (block.nd == 2)
      return (block.ni+1)*(block.nj+1);
    else
      return (block.ni+1)*(block.nj+1)*(block.nk+1);
  }
};

////////////////////////////////////////////////////////////////////////////////
/// Class to iterate over vertices
////////////////////////////////////////////////////////////////////////////////
struct vertex_iterator_t {
  const logical_vertex_t * vertices = nullptr;

  int_t i=0, j=0, k=0;
  
  vertex_iterator_t() = default;

  vertex_iterator_t(const logical_vertex_t & vs) : vertices(&vs) 
  {}

  vertex_iterator_t(const logical_vertex_t & vs, int_t ii, int_t jj, int_t kk) :
    vertices(&vs), i(ii), j(jj), k(kk) {}

  void reset() {
    i = 0;
    j = 0;
    k = 0;
  }
  
  bool next() {
    ++i;
    if (i >= vertices->ni) {
      i = 0;
      ++j;
      if (j >= vertices->nj) {
        j = 0;
        ++k;
        if (k >= vertices->nk) {
          return false;
        }
      }
    }
    return true;
  }
  
  bool operator<(const vertex_iterator_t & other)
  { 
    if (i != other.i) return i<other.i;
    if (j != other.j) return j<other.j;
    return k < other.k;
  }

  int_t id() const { return vertices->id(i, j, k); }

  void pos(int_t * vals) {
    vals[0] = i;
    vals[1] = j;
    vals[2] = k;
  }

  vertex_iterator_t offset(int_t ioff, int_t joff, int_t koff)
  { return vertex_iterator_t(*vertices, i+ioff, j+joff, k+koff); }
  
  int_t offset_id(int_t ioff, int_t joff, int_t koff)
  { return vertices->id(i+ioff, j+joff, k+koff); }
  
  
  vertex_iterator_t offset(int_t off, int_t dir)
  {
    switch (dir) {
      case 0:   return vertex_iterator_t(*vertices, i+off, j, k);
      case 1:   return vertex_iterator_t(*vertices, i, j+off, k);
      default:  return vertex_iterator_t(*vertices, i, j, k+off);
    }
  }
  
  int_t offset_id(int_t off, int_t dir)
  {
    switch (dir) {
      case 0:   return vertices->id(i+off, j, k);
      case 1:   return vertices->id(i, j+off, k);
      default:  return vertices->id(i, j, k+off);
    }
  }


  bool is_iboundary(bool is_hi) const
  { return is_hi ? i==vertices->ni-1 : i==0; }
  bool is_iboundary() const
  { return i==vertices->ni-1 || i==0; }

  bool is_jboundary(bool is_hi) const
  { return is_hi ? j==vertices->nj-1 : j==0; }
  bool is_jboundary() const
  { return j==vertices->nj-1 || j==0; }

  bool is_kboundary(bool is_hi) const
  { return is_hi ? k==vertices->nk-1 : k==0; }
  bool is_kboundary() const
  { return k==vertices->nk-1 || k==0; }

  bool is_boundary(bool is_hi, int_t dir) const {
    switch (dir) {
      case 0:   return is_hi ? i==vertices->ni-1 : i==0;
      case 1:   return is_hi ? j==vertices->nj-1 : j==0;
      default:  return is_hi ? k==vertices->nk-1 : k==0;
    }
  }
  
  bool is_boundary(int_t dir) const {
    switch (dir) {
      case 0:   return i==vertices->ni-1 || i==0;
      case 1:   return j==vertices->nj-1 || j==0;
      default:  return k==vertices->nk-1 || k==0;
    }
  }
  
  bool is_boundary() const {
    auto nd = vertices->block.nd;
    return 
      (nd>0 && (i==vertices->ni-1 || i==0)) ||
      (nd>1 && (j==vertices->nj-1 || j==0)) ||
      (nd>2 && (k==vertices->nk-1 || k==0));
  }
  
  bool is_sector_boundary(const int_t * sector)
  {
    bool is_sector_bnd = true;
    auto nd = vertices->block.nd;
    for (int_t dim=0; dim<nd; ++dim) {
      if (sector[dim]) {
        auto is_bnd = sector[dim]<0 ? is_boundary(0,dim) : is_boundary(1,dim);
        is_sector_bnd = is_sector_bnd && is_bnd;
      }
    } // dim
    return is_sector_bnd;
  }
};
  

////////////////////////////////////////////////////////////////////////////////
/// Class to iterate over edges
////////////////////////////////////////////////////////////////////////////////
struct logical_edge_t {
  const logical_block_t & block;
  
  int_t dim=0;
  int_t start = 0, end = 0;
  int_t strides[3] = {0, 0, 0};
  int_t ni=0, nj=0, nk=0;

  logical_edge_t(const logical_block_t & blk, int_t d) :
    block(blk), dim(d),
    ni(block.ni), nj(block.nj), nk(block.nk)
  {
    auto nd = blk.nd;
    
    auto nip = block.ni+1;
    auto njp = block.nj+(nd>1);
    auto nkp = block.nk+(nd>2);

    switch (dim) {
      case 0:
        strides[0] = 1;
        strides[1] = ni;
        strides[2] = ni*njp;
        end = ni*njp*nkp;
        nj = njp;
        nk = nkp;
        break;
      case 1:
        strides[0] = 1;
        strides[1] = nip;
        strides[2] = nip*nj;
        start = ni*njp*nkp;
        end = start + nip*nj*nkp;
        ni = nip;
        nk = nkp;
        break;
      default:
        strides[0] = 1;
        strides[1] = nip;
        strides[2] = nip*njp;
        start = ni*njp*nkp + nip*nj*nkp;
        end = start + nip*njp*nk;
        ni = nip;
        nj = njp;
    };
  }

  int_t id(int_t i, int_t j, int_t k) const
  { return start + i*strides[0] + j*strides[1] + k*strides[2]; }
  
  edge_iterator_t begin() const;
  int size() const { return end - start; }

};

////////////////////////////////////////////////////////////////////////////////
/// Underying edge iterator
////////////////////////////////////////////////////////////////////////////////
struct edge_iterator_t {

  const logical_edge_t * edges = nullptr;

  int_t i=0, j=0, k=0;

  edge_iterator_t() = default;
  
  edge_iterator_t(const logical_edge_t & es) : edges(&es) {}

  edge_iterator_t(const logical_edge_t & es, int_t ii, int_t jj, int_t kk) :
    edges(&es), i(ii), j(jj), k(kk)
  {}

  void reset() {
    i = 0;
    j = 0;
    k = 0;
  }
  
  bool next() {
    ++i;
    if (i >= edges->ni) {
      i = 0;
      ++j;
      if (j >= edges->nj) {
        j = 0;
        ++k;
        if (k >= edges->nk) {
          return false;
        }
      }
    }
    return true;
  }

  edge_iterator_t & operator++() {
    next();
    return *this;
  }

  bool operator<(const edge_iterator_t & other)
  { 
    if (i != other.i) return i<other.i;
    if (j != other.j) return j<other.j;
    return k < other.k;
  }

  int_t id() const { return edges->id(i, j, k); }

  void vertex_ids(const logical_vertex_t &vertices, int_t * vs) const
  {
    vs[0] = vertices.id(i, j, k);
    vs[1] = vertices.id(i+(edges->dim==0), j+(edges->dim==1), k+(edges->dim==2));
  }

};
  

////////////////////////////////////////////////////////////////////////////////
/// Class to iterate over faces
////////////////////////////////////////////////////////////////////////////////
struct logical_face_t {
  const logical_block_t & block;
  
  int_t dim=0;
  int_t start = 0, end = 0;
  int_t strides[3] = {0, 0, 0};

  logical_face_t(const logical_block_t & blk, int_t d) :
    block(blk), dim(d)
  {
    switch (dim) {
      case 0:
        strides[0] = 1;
        strides[1] = block.ni+1;
        strides[2] = (block.ni+1)*block.nj;
        end = (block.ni+1)*block.nj*block.nk;
        break;
      case 1:
        strides[0] = 1;
        strides[1] = block.ni;
        strides[2] = block.ni*(block.nj+1);
        start = (block.ni+1)*block.nj*block.nk;
        end = start + block.ni*(block.nj+1)*block.nk;
        break;
      default:
        strides[0] = 1;
        strides[1] = block.ni;
        strides[2] = block.ni*block.nj;
        start = (block.ni+1)*block.nj*block.nk + block.ni*(block.nj+1)*block.nk;
        end = start + block.ni*block.nj*(block.nk+1);
    };
  }

  int_t id(int_t i, int_t j, int_t k) const
  { return start + i*strides[0] + j*strides[1] + k*strides[2]; }
  
  face_iterator_t begin() const;
  int size() const { return end - start; }

  face_iterator_t at(int_t i, int_t j, int_t k) const;

  face_iterator_t at(int_t id) const;

};

////////////////////////////////////////////////////////////////////////////////
/// Underying face iterator
////////////////////////////////////////////////////////////////////////////////
struct face_iterator_t {

  const logical_face_t * faces = nullptr;

  int_t i=0, j=0, k=0;

  face_iterator_t() = default;

  face_iterator_t(const logical_face_t & fs) : faces(&fs) {}

  face_iterator_t(const logical_face_t & fs, int_t ii, int_t jj, int_t kk) :
    faces(&fs), i(ii), j(jj), k(kk)
  {}

  int_t dim() const { return faces->dim; }

  void reset() {
    i = 0;
    j = 0;
    k = 0;
  }
  
  bool next() {

    ++i;
    if (i >= faces->block.ni + (faces->dim==0)) {
      i = 0;
      ++j;
      if (j >= faces->block.nj + (faces->dim==1)) {
        j = 0;
        ++k;
        if (k >= faces->block.nk + (faces->dim==2)) {
          return false;
        }
      }
    }
    return true;

  }

  face_iterator_t & operator++() {
    next();
    return *this;
  }

  bool operator<(const face_iterator_t & other)
  { 
    if (i != other.i) return i<other.i;
    if (j != other.j) return j<other.j;
    return k < other.k;
  }

  int_t id() const { return faces->id(i, j, k); }

  cell_iterator_t left_cell(const logical_cell_t &) const;
  cell_iterator_t right_cell(const logical_cell_t &) const;
  
  int_t left_cell_id(const logical_cell_t & cells) const;
  int_t right_cell_id(const logical_cell_t & cells) const;
  
  void vertex_ids(const logical_vertex_t &vertices, int_t * vs) const
  {
    if (faces->block.nd == 1) {
      vs[0] = vertices.id(i, j, k);
    }
    else if (faces->block.nd == 2) {
      vs[0] = vertices.id(i, j,   k);
      vs[1] = vertices.id(i+(faces->dim==1), j+(faces->dim==0), k);
    }
    else if (faces->block.nd == 3) {
      bool ifa = faces->dim==0;
      bool jfa = faces->dim==1;
      bool kfa = faces->dim==2;
      vs[0] = vertices.id(i, j, k);
      vs[1] = vertices.id(i+(jfa||kfa), j, k+(ifa));
      vs[2] = vertices.id(i+(jfa||kfa), j+(ifa||kfa), k+(ifa||jfa));
      vs[3] = vertices.id(i, j+(ifa||kfa), k+(jfa));
    }
  }
  
  void edge_ids(const logical_edge_t &edges, int_t * es) const
  {
    if (faces->block.nd == 3 && edges.dim != faces->dim) {
      auto ei = edges.dim==0;
      auto ej = edges.dim==1;
      auto ek = edges.dim==2;
      auto fi = faces->dim==0;
      auto fj = faces->dim==1;
      auto fk = faces->dim==2;
      es[0] = edges.id( i, j, k );
      es[1] = edges.id( 
          i + (fj&&ek) + (fk&&ej),
          j + (fi&&ek) + (fk&&ei),
          k + (fi&&ej) + (fj&&ei) );
    }
  }
  
  void edges(const logical_edge_t &edges, edge_iterator_t * es) const
  {
    if (faces->block.nd == 3 && edges.dim != faces->dim) {
      auto ei = edges.dim==0;
      auto ej = edges.dim==1;
      auto ek = edges.dim==2;
      auto fi = faces->dim==0;
      auto fj = faces->dim==1;
      auto fk = faces->dim==2;
      es[0] = edge_iterator_t( edges, i, j, k );
      es[1] = edge_iterator_t( 
          edges,
          i + (fj&&ek) + (fk&&ej),
          j + (fi&&ek) + (fk&&ei),
          k + (fi&&ej) + (fj&&ei) );
    }
  }
  
  
  bool is_boundary() 
  {
    switch (faces->dim) {
      case(0): return i==0 || i==faces->block.ni;
      case(1): return j==0 || j==faces->block.nj;
      default: return k==0 || k==faces->block.nk;
    }
  }
  bool is_boundary(bool is_hi) 
  {
    switch (faces->dim) {
      case(0): return is_hi ? i==faces->block.ni : i==0;
      case(1): return is_hi ? j==faces->block.nj : j==0;
      default: return is_hi ? k==faces->block.nk : k==0;
    }
  }
  
  bool is_iboundary() 
  { return faces->dim==0 && (i==0 || i==faces->block.ni); }
  bool is_iboundary(bool is_hi) 
  { return is_hi ? faces->dim==0&&i==faces->block.ni : faces->dim==0&&i==0; }
  
  bool is_jboundary() 
  { return faces->dim==1 && (j==0 || j==faces->block.nj); }
  bool is_jboundary(bool is_hi) 
  { return is_hi ? faces->dim==1&&j==faces->block.nj : faces->dim==1&&j==0; }
  
  bool is_kboundary() 
  { return faces->dim==2 && (k==0 || k==faces->block.nk); }
  bool is_kboundary(bool is_hi) 
  { return is_hi ? faces->dim==2&&k==faces->block.nk : faces->dim==2&&k==0; }

};
  
////////////////////////////////////////////////////////////////////////////////
/// Class to iterate over cells
////////////////////////////////////////////////////////////////////////////////
struct logical_cell_t {
 
  const logical_block_t & block;

  int_t strides[3] = {0, 0, 0};
 
  logical_cell_t(const logical_block_t & blk) :
    block(blk)
  {
    strides[0] = 1;
    strides[1] = block.ni;
    strides[2] = block.ni * block.nj;
  }

  int_t id(int_t i, int_t j, int_t k) const
  { return i*strides[0] + j*strides[1] + k*strides[2]; }
  
  cell_iterator_t at(int_t pos) const;
  
  cell_iterator_t at(int_t i, int_t j, int_t k) const;

  cell_iterator_t begin() const;
  
  int_t size() const
  { return block.ni*block.nj*block.nk;  }

};
  
////////////////////////////////////////////////////////////////////////////////
/// Underying cell iterator
////////////////////////////////////////////////////////////////////////////////
struct cell_iterator_t {

  const  logical_cell_t & cells;
  int_t i=0, j=0, k=0;

  cell_iterator_t(const logical_cell_t & cs) : cells(cs) {}

  cell_iterator_t(const logical_cell_t & cs, int_t ii, int_t jj, int_t kk) :
    cells(cs), i(ii), j(jj), k(kk) {}

  void reset() {
    i = 0;
    j = 0;
    k = 0;
  }

  bool next() {
    ++i;
    if (i >= cells.block.ni) {
      i = 0;
      ++j;
      if (j >= cells.block.nj) {
        j = 0;
        ++k;
        if (k >= cells.block.nk) {
          return false;
        }
      }
    }
    return true;
  }
  
  bool next(int_t idir, int_t jdir, int_t kdir) {
    i+=idir;
    if ( (idir>0 && i>=cells.block.ni) || (idir<0 && i<0) ) {
      i = idir>0 ? 0 : cells.block.ni-1;
      j+=jdir;
      if ( (jdir>0 && j>=cells.block.nj) || (jdir<0 && j<0) ) {
        j = jdir>0 ? 0 : cells.block.nj-1;
        k+=kdir;
        if ( (kdir>0 && k>=cells.block.nk) || (kdir<0 && k<0) ) {
          return false;
        }
      }
    }
    return true;
  }
  
  bool operator<(const cell_iterator_t & other)
  { 
    if (i != other.i) return i<other.i;
    if (j != other.j) return j<other.j;
    return k < other.k;
  }

  int_t id() const { return cells.id(i, j, k); }

  bool is_iboundary(bool is_hi) const
  { return is_hi ? i==cells.block.ni-1 : i==0; }
  bool is_iboundary() const
  { return i==cells.block.ni-1 || i==0; }

  bool is_jboundary(bool is_hi) const
  { return is_hi ? j==cells.block.nj-1 : j==0; }
  bool is_jboundary() const
  { return j==cells.block.nj-1 || j==0; }

  bool is_kboundary(bool is_hi) const
  { return is_hi ? k==cells.block.nk-1 : k==0; }
  bool is_kboundary() const
  { return k==cells.block.nk-1 || k==0; }

  bool is_boundary(bool is_hi, int_t dir) const {
    switch (dir) {
      case 0:   return is_hi ? i==cells.block.ni-1 : i==0;
      case 1:   return is_hi ? j==cells.block.nj-1 : j==0;
      default:  return is_hi ? k==cells.block.nk-1 : k==0;
    }
  }
  
  bool is_boundary(int_t dir) const {
    switch (dir) {
      case 0:   return i==cells.block.ni-1 || i==0;
      case 1:   return j==cells.block.nj-1 || j==0;
      default:  return k==cells.block.nk-1 || k==0;
    }
  }
  
  bool is_boundary() const {
    auto nd = cells.block.nd;
    return 
      (nd>0 && (i==cells.block.ni-1 || i==0)) ||
      (nd>1 && (j==cells.block.nj-1 || j==0)) ||
      (nd>2 && (k==cells.block.nk-1 || k==0));
  }

  cell_iterator_t offset(int_t ioff, int_t joff, int_t koff) const
  { return cell_iterator_t(cells, i+ioff, j+joff, k+koff); }
  
  int_t offset_id(int_t ioff, int_t joff, int_t koff) const
  { return cells.id(i+ioff, j+joff, k+koff); }

  cell_iterator_t offset(int_t off, int_t dir) const
  {
    switch (dir) {
      case 0:   return cell_iterator_t(cells, i+off, j, k);
      case 1:   return cell_iterator_t(cells, i, j+off, k);
      default:  return cell_iterator_t(cells, i, j, k+off);
    }
  }
  
  int_t offset_id(int_t off, int_t dir) const
  {
    switch (dir) {
      case 0:   return cells.id(i+off, j, k);
      case 1:   return cells.id(i, j+off, k);
      default:  return cells.id(i, j, k+off);
    }
  }

  face_iterator_t face(const logical_face_t &faces, int_t off) const
  {
    switch (faces.dim) {
      case 0:   return face_iterator_t(faces, i+off, j, k);
      case 1:   return face_iterator_t(faces, i, j+off, k);
      default:  return face_iterator_t(faces, i, j, k+off);
    }
  }
  
  int_t face_id(const logical_face_t &faces, int_t off) const
  {
    switch (faces.dim) {
      case 0:   return faces.id(i+off, j, k);
      case 1:   return faces.id(i, j+off, k);
      default:  return faces.id(i, j, k+off);
    }
  }

  vertex_iterator_t vertex(
      const logical_vertex_t &vertices,
      int_t ioff,
      int_t joff,
      int_t koff) const
  { return vertex_iterator_t(vertices, i+ioff, j+joff, k+koff); }
  
  int_t vertex_id(
      const logical_vertex_t &vertices,
      int_t ioff,
      int_t joff,
      int_t koff) const
  { return vertices.id(i+ioff, j+joff, k+koff); }
  
  void vertices(const logical_vertex_t &vertices, vertex_iterator_t * vs) const
  {
    if (cells.block.nd == 1) {
      vs[0] = vertex_iterator_t(vertices, i,   j, k);
      vs[1] = vertex_iterator_t(vertices, i+1, j, k);
    }
    else if (cells.block.nd == 2) {
      vs[0] = vertex_iterator_t(vertices, i,   j,   k);
      vs[1] = vertex_iterator_t(vertices, i+1, j,   k);
      vs[2] = vertex_iterator_t(vertices, i+1, j+1, k);
      vs[3] = vertex_iterator_t(vertices, i,   j+1, k);
    }
    else if (cells.block.nd == 3) {
      vs[0] = vertex_iterator_t(vertices, i,   j,   k);
      vs[1] = vertex_iterator_t(vertices, i+1, j,   k);
      vs[2] = vertex_iterator_t(vertices, i+1, j+1, k);
      vs[3] = vertex_iterator_t(vertices, i,   j+1, k);
      vs[4] = vertex_iterator_t(vertices, i,   j,   k+1);
      vs[5] = vertex_iterator_t(vertices, i+1, j,   k+1);
      vs[6] = vertex_iterator_t(vertices, i+1, j+1, k+1);
      vs[7] = vertex_iterator_t(vertices, i,   j+1, k+1);
    }
  }
  
  void vertex_ids(const logical_vertex_t &vertices, int_t * vs) const
  {
    if (cells.block.nd == 1) {
      vs[0] = vertices.id(i, j, k);
      vs[1] = vertices.id(i+1, j, k);
    }
    else if (cells.block.nd == 2) {
      vs[0] = vertices.id(i,   j,   k);
      vs[1] = vertices.id(i+1, j,   k);
      vs[2] = vertices.id(i+1, j+1, k);
      vs[3] = vertices.id(i,   j+1, k);
    }
    else if (cells.block.nd == 3) {
      vs[0] = vertices.id(i,   j,   k);
      vs[1] = vertices.id(i+1, j,   k);
      vs[2] = vertices.id(i+1, j+1, k);
      vs[3] = vertices.id(i,   j+1, k);
      vs[4] = vertices.id(i,   j,   k+1);
      vs[5] = vertices.id(i+1, j,   k+1);
      vs[6] = vertices.id(i+1, j+1, k+1);
      vs[7] = vertices.id(i,   j+1, k+1);
    }
  }

  void edges(const logical_edge_t &edges, edge_iterator_t * es) const
  {
    if (cells.block.nd == 1) {
      es[0] = edge_iterator_t(edges, i, j, k);
    }
    else if (cells.block.nd == 2) {
      auto isi = edges.dim==0;
      auto isj = edges.dim==1;
      es[0] = edge_iterator_t(edges, i,     j,       k);
      es[1] = edge_iterator_t(edges, i+isj, j+isi,   k);
    }
    else if (cells.block.nd == 3) {
      auto isi = edges.dim==0;
      auto isj = edges.dim==1;
      auto isk = edges.dim==2;
      es[0] = edge_iterator_t(edges, i, j, k);
      es[1] = edge_iterator_t(edges, i+(isj||isk), j+isi, k);
      es[2] = edge_iterator_t(edges, i+(isj||isk), j+(isi||isk), k+(isi||isj));
      es[3] = edge_iterator_t(edges, i, j+isk, k+(isi||isj));
    }
  }

  void edge_ids(const logical_edge_t &edges, int_t * es) const
  {
    if (cells.block.nd == 1) {
      es[0] = edges.id(i, j, k);
    }
    else if (cells.block.nd == 2) {
      auto isi = edges.dim==0;
      auto isj = edges.dim==1;
      es[0] = edges.id( i,     j,       k);
      es[1] = edges.id( i+isj, j+isi,   k);
    }
    else if (cells.block.nd == 3) {
      auto isi = edges.dim==0;
      auto isj = edges.dim==1;
      auto isk = edges.dim==2;
      es[0] = edges.id( i, j, k);
      es[1] = edges.id( i+(isj||isk), j+isi, k);
      es[2] = edges.id( i+(isj||isk), j+(isi||isk), k+(isi||isj));
      es[3] = edges.id( i, j+isk, k+(isi||isj));
    }
  }

}; // iterator

////////////////////////////////////////////////////////////////////////////////
/// vertex declarations
////////////////////////////////////////////////////////////////////////////////
inline vertex_iterator_t logical_vertex_t::begin() const
{ return vertex_iterator_t(*this); }

  
////////////////////////////////////////////////////////////////////////////////
/// Face declarations
////////////////////////////////////////////////////////////////////////////////
inline
cell_iterator_t
face_iterator_t::left_cell(const logical_cell_t & cells) const
{
  switch (faces->dim) {
    case (0): return cell_iterator_t(cells, i-1, j, k);
    case (1): return cell_iterator_t(cells, i, j-1, k);
    default:  return cell_iterator_t(cells, i, j, k-1);
  }
}
  
inline
int_t face_iterator_t::left_cell_id(const logical_cell_t & cells) const
{
  switch (faces->dim) {
    case (0): return cells.id(i-1, j, k);
    case (1): return cells.id(i, j-1, k);
    default:  return cells.id(i, j, k-1);
  }
}

inline face_iterator_t logical_face_t::begin() const
{ return face_iterator_t(*this); }

inline
cell_iterator_t
face_iterator_t::right_cell(const logical_cell_t & cells) const
{ return cell_iterator_t(cells, i, j, k); }

inline
int_t face_iterator_t::right_cell_id(const logical_cell_t & cells) const
{ return cells.id(i, j, k); }

inline face_iterator_t logical_face_t::at(int_t i, int_t j, int_t k) const
{ return face_iterator_t(*this, i, j, k); }

inline face_iterator_t logical_face_t::at(int_t id) const
{ 
  auto pos = id - start;

  if (pos<0) return face_iterator_t();

  switch (block.nd) {
    case (0):
      return face_iterator_t(*this, pos, 0, 0);
    case (1): {
      int_t j = pos / strides[1];
      int_t i = pos % strides[1]; 
      return face_iterator_t(*this, i, j, 0);
    }
    default: {
      int_t k = pos / strides[2];
      pos -= k*strides[2];
      int_t j = pos / strides[1];
      int_t i = pos % strides[1]; 
      return face_iterator_t(*this, i, j, k);
    }
  }
}

inline face_iterator_t at(const logical_face_t * faces, int_t id)
{
  auto nd = faces->block.nd;
  for (int_t d=0; d<nd; ++d) {
    if (id < faces[d].end)
      return faces[d].at(id);
  }
  return {};
}

////////////////////////////////////////////////////////////////////////////////
/// cell declarations
////////////////////////////////////////////////////////////////////////////////
inline cell_iterator_t logical_cell_t::begin() const
{ return cell_iterator_t(*this); }


inline cell_iterator_t logical_cell_t::at(int_t pos) const
{
  switch (block.nd) {
    case (0):
      return cell_iterator_t(*this, pos, 0, 0);
    case (1): {
      int_t j = pos / strides[1];
      int_t i = pos % strides[1]; 
      return cell_iterator_t(*this, i, j, 0);
    }
    default: {
      int_t k = pos / strides[2];
      pos -= k*strides[2];
      int_t j = pos / strides[1];
      int_t i = pos % strides[1]; 
      return cell_iterator_t(*this, i, j, k);
    }
  };
}
  

inline cell_iterator_t logical_cell_t::at(int_t i, int_t j, int_t k) const
{ return cell_iterator_t(*this, i, j, k); }


////////////////////////////////////////////////////////////////////////////////
/// Block declarations
////////////////////////////////////////////////////////////////////////////////
inline logical_vertex_t logical_block_t::vertices() const
{ return logical_vertex_t(*this); }

inline logical_edge_t logical_block_t::edges(int_t dim) const
{ return logical_edge_t(*this, dim); }
  
inline logical_face_t logical_block_t::faces(int_t dim) const
{ return logical_face_t(*this, dim); }

inline logical_cell_t logical_block_t::cells() const
{ return logical_cell_t(*this); }

} // namespace

#endif
