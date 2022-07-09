#ifndef CONNECT_HPP
#define CONNECT_HPP

#include "utils/crs.hpp"

namespace prl {

//==============================================================================
//! \brief Create connectivity among cells, edges, and vertices from known
//! cell-to-vertex connectivity.
//==============================================================================
template<typename T, typename FUNCTION>
void
build(
    const crs_t<T,int_t> & cell_to_vertex,
    crs_t<T,int_t> & cell_to_edge,
    crs_t<T,int_t> & edge_to_vertex,
    std::map<std::vector<T>, T> & sorted_vertices_to_edges,
    FUNCTION && build_edges_from_vertices)
{

  // some temporary storage
  crs_t<T, int_t> new_edges;
  std::vector<T> sorted_vs;
  std::vector<T> these_edges;
  
  // loop over cells, adding all of their edges to the table
  for (const auto & these_verts : cell_to_vertex) {

    // clear this cells edges
    these_edges.clear();

    // now build the edges for the cell
    new_edges.clear();
    auto make_new = std::forward<FUNCTION>(build_edges_from_vertices)(these_verts, new_edges);

    // now look for exsiting vertex pairs in the edge-to-vertex master list
    // or add a new edge to the list.  add the matched edge id to the
    // cell-to-edge list
    for (const auto & vs : new_edges) {
      // sort the vertices
      sorted_vs.assign(vs.begin(), vs.end());
      std::sort(sorted_vs.begin(), sorted_vs.end());
      // if we dont find the edge
      auto it = sorted_vertices_to_edges.find(sorted_vs);
      if (it == sorted_vertices_to_edges.end())
      {
        if (make_new) {
          // add to the local reverse map
          auto edgeid = sorted_vertices_to_edges.size();
          sorted_vertices_to_edges.emplace(sorted_vs, edgeid);
          // add to the original sorted and unsorted maps
          edge_to_vertex.push_back(vs);
          // add to the list of edges
          these_edges.push_back(edgeid);
        } // make
      }
      else {
        // if we do find the edge
        // just add the id to the list of edges
        these_edges.push_back(it->second);
      }
    }

    // now add this cells edges
    cell_to_edge.push_back( these_edges );

  }  // for
}  // build_connectivity

//==============================================================================
//! \brief Transpose a connectivity array.
//==============================================================================
template<typename T>
void
transpose(const crs_t<T, int_t> & in, crs_t<T, int_t> & out) {
  
  // determine offsets
  size_t num_from = in.size();

  const auto & in_offsets = in.offsets;
  const auto & in_indices = in.indices;
  size_t num_to = *std::max_element(in_indices.begin(), in_indices.end()) + 1;
  
  auto & out_offsets = out.offsets;
  out_offsets.clear();
  out_offsets.resize(num_to+1, 0);
  
  for (size_t from = 0; from < num_from; ++from) {
    for (auto i=in_offsets[from]; i<in_offsets[from+1]; ++i) {
      auto to = in_indices[i];
      out_offsets[to+1]++;
    }
  }

  out_offsets[0] = 0;
  for (size_t to = 0; to < num_to; ++to)
    out_offsets[to+1] += out_offsets[to];


  // determine final connectivity
  auto & out_indices = out.indices;
  out_indices.clear();
  out_indices.resize( out_offsets.back() );

  std::vector<int> pos(num_to, 0);

  for (size_t from = 0; from < num_from; ++from) {
    for (auto i=in_offsets[from]; i<in_offsets[from+1]; ++i) {
      auto to = in_indices[i];
      out_indices[ out_offsets[to] + pos[to] ] = from;
      pos[to]++;
    }
  }
}

//==============================================================================
//! \brief Intersect two connectivity arrays; that is, given X-to-Y and Y-to-Z
//! connectivity, build X-to-Z connectivity.
//==============================================================================
template<typename T>
void
intersect(
    const crs_t<T, int_t> & cell_to_face,
    const crs_t<T, int_t> & face_to_edge,
    crs_t<T, int_t> & cell_to_edge)
{

  // some temporary storage
  std::vector<T> these_edges;

  // loop over cells, collecting their face edges
  for (const auto & these_faces : cell_to_face) {
    // reference the storage for this cells edges
    these_edges.clear();

    // now add the edges of this cell
    for (auto f : these_faces)
      for (auto e : face_to_edge.at(f)) {
        // set insertion is not the quickest, but maintains order
        auto it = std::find( these_edges.begin(), these_edges.end(), e );
        if ( it == these_edges.end() ) these_edges.push_back(e);
      }

    // sort them and remove the non-unique ones
    // ... This mightbe quicker but does not maintain order
    // std::sort(these_edges.begin(), these_edges.end());
    // auto last = std::unique(these_edges.begin(), these_edges.end());
    // these_edges.erase(last, these_edges.end());

    // now add this cells edges
    cell_to_edge.push_back( these_edges );

  }  // cells

}  // intersect

//==============================================================================
//! \brief Determine neighbors
//==============================================================================
template<typename T>
void
neighbors(
    const crs_t<T, int_t> & conn,
    const crs_t<T, int_t> & rconn,
    crs_t<T, int_t> & neigh)
{

  // some temporary storage
  std::vector<T> ents;

  auto & offsets = neigh.offsets;
  auto & indices = neigh.indices;

  offsets.clear();
  indices.clear();

  offsets.emplace_back(0);

  // loop over cells, collecting their face edges
  T cnt = 0;
  for (const auto & to_ents : conn) {
    // reference the storage for this cells edges
    ents.clear();

    // now add the edges of this cell
    for (auto to : to_ents)
      for (auto from : rconn.at(to))
        if (cnt != from)
          ents.push_back(from);

    // sort them and remove the non-unique ones
    // ... This mightbe quicker but does not maintain order
    std::sort(ents.begin(), ents.end());
    auto last = std::unique(ents.begin(), ents.end());
    
    indices.insert(indices.end(), ents.begin(), last);
    offsets.emplace_back(indices.size());

    cnt++;

  }  // cells

}  // intersect

} // namespace

#endif

