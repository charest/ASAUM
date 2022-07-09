#ifndef EXO_WRITER_HPP
#define EXO_WRITER_HPP

#include "exo_types.hpp"
#include "shape.hpp"
#include "utils/crs.hpp"
#include "utils/errors.hpp"

#include <exodusII.h>

#include <algorithm>
#include <cstring>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace prl {

class lua_result_t;
class mpi_comm_t;
struct outvar_t;

////////////////////////////////////////////////////////////////////////////////
/// an exodus reader class
////////////////////////////////////////////////////////////////////////////////
class exo_writer_t {
public:
  
  /// the precision type
  enum class precision_t { sp, dp, def };
  
  /// the file id
  int exoid_ = -1;

  /// number of dimensions
  int num_dims_ = 0;

  /// the precision
  precision_t precision_ = precision_t::def;

  // variables
  struct varinfo_t {
    int_t id;
    int_t components;
  };
  std::map<std::string, varinfo_t> var_map_;

public:
  
  /// constructors
  exo_writer_t() = default;
  exo_writer_t(lua_result_t input, mpi_comm_t & comm);

  int open(const std::string & name);
  int close();

  bool is_int64();

  // write parameters to exodus
  static ex_init_params make_params(const char * label = nullptr);
  int write_params(const ex_init_params & exo_params);
  
  // write fields
  int write_field_names(const std::vector<outvar_t> & vars);
  int write_field(const real_t * data, size_t n, const std::string & name, int_t blk);

  int write_time(size_t time_step, real_t time);
  
  int write_side_set_param(int_t ss_id, size_t num_sides);
  int write_side_set_names(const std::vector<std::string> & names);


  //============================================================================
  //! \brief get the shape string
  //============================================================================
  static const char * geom_type(shape_t shape) {
    switch (shape) {
      case shape_t::line:    return "bar";
      case shape_t::tri:     return "tri3";
      case shape_t::quad:    return "quad4";
      case shape_t::polygon: return "nsided";
      case shape_t::tet:     return "tet4";
      case shape_t::hex:     return "hex8";
      case shape_t::polyhedron: return "nfaced";
      default:
        THROW_ERROR("Unknown shape type: " << static_cast<int>(shape));
    }
    return nullptr;
  }
  
  
  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename ENTITY_CONN>
  int write_block(
      size_t blk_id,
      const std::string & name,
      ex_entity_type entity_type,
      const std::string & entity_description,
      size_t num_entities,
      ENTITY_CONN && entity_conn)
  {

    // some type aliases
    using ex_index_t = U;

    //-------------------------------------------------------------------------
    if ( entity_description == "nfaced" || entity_description == "nsided" )
    {
      // check if face data is provided instead of node data
      auto is_face_data = (num_dims_ == 3 && entity_type == EX_ELEM_BLOCK);

      // guess how many elements are in each connectivity slot
      auto num_conn_guess = num_dims_ * num_dims_;

      // build the connectivitiy list for the block
      std::vector<ex_index_t> entity_nodes;
      std::vector<int> entity_node_counts;
      entity_nodes.reserve(num_entities * num_conn_guess);
      entity_node_counts.reserve(num_entities);

      // temporary storage for each connecitivity slot
      std::vector<ex_index_t> temp_entity_nodes;
      temp_entity_nodes.reserve(num_conn_guess);

      for (size_t e = 0; e < num_entities; ++e) {
        // get the connectivity from the user-provided function
        std::forward<ENTITY_CONN>(entity_conn)(e, temp_entity_nodes);
        // store them in the master list ( exodus wants 1-index arrays )
        for (auto i : temp_entity_nodes)
          entity_nodes.emplace_back(i + 1);
        entity_node_counts.push_back(temp_entity_nodes.size());
        // reset temp storage
        temp_entity_nodes.clear();
      }

      // the total size needed to hold the element connectivity
      ex_index_t num_nodes_this_blk = entity_nodes.size();
      ex_index_t num_entries_this_blk = entity_node_counts.size();

      // set the block header
      ex_index_t num_attr_per_entry = 0;
      ex_index_t num_nodes_per_entry = is_face_data ? 0 : num_nodes_this_blk;
      ex_index_t num_edges_per_entry = 0;
      ex_index_t num_faces_per_entry = is_face_data ? num_nodes_this_blk : 0;
      auto status = ex_put_block(
          exoid_, entity_type, blk_id, entity_description.c_str(), num_entries_this_blk,
          num_nodes_per_entry, num_edges_per_entry, num_faces_per_entry,
          num_attr_per_entry);
      if (status) return status;

      // write the block name
      status = ex_put_name(exoid_, entity_type, blk_id, name.c_str());
      if (status) return status;

      // write connectivity
      auto node_conn = is_face_data ? nullptr : entity_nodes.data();
      auto elem_edge_conn = nullptr;
      auto elem_face_conn = is_face_data ? entity_nodes.data() : nullptr;
      status = ex_put_conn(
          exoid_, entity_type, blk_id, node_conn, elem_edge_conn, elem_face_conn);
      if (status) return status;

      // write counts
      status = ex_put_entity_count_per_polyhedra(
          exoid_, entity_type, blk_id, entity_node_counts.data());
      if (status) return status;

    }

    //-------------------------------------------------------------------------
    else {
      
      std::vector<ex_index_t> entity_nodes, temp_entity_nodes;
      
      for (size_t e = 0; e < num_entities; ++e) {
        temp_entity_nodes.clear();
        // get the connectivity from the user-provided function
        std::forward<ENTITY_CONN>(entity_conn)(e, temp_entity_nodes);
        // store them in the master list ( exodus wants 1-index arrays )
        for (auto i : temp_entity_nodes)
          entity_nodes.emplace_back(i + 1);
      }

      // write the block header
      ex_index_t num_entries_this_blk = num_entities;
      ex_index_t num_nodes_per_entry = temp_entity_nodes.size();
      ex_index_t num_edges_per_entry = 0;
      ex_index_t num_faces_per_entry = 0;
      ex_index_t num_attr_per_entry = 0;
      
      auto status = ex_put_block(
          exoid_, entity_type, blk_id, entity_description.c_str(),
          num_entries_this_blk, num_nodes_per_entry, num_edges_per_entry,
          num_faces_per_entry, num_attr_per_entry);
      if (status) return status;

      // write the block name
      status = ex_put_name(exoid_, entity_type, blk_id, name.c_str());
      if (status) return status;

      // write connectivity
      auto node_conn = entity_nodes.data();
      auto elem_edge_conn = nullptr;
      auto elem_face_conn = nullptr;
      status = ex_put_conn(
          exoid_, entity_type, blk_id, node_conn, elem_edge_conn, elem_face_conn);
      if (status) return status;
    } // type

    return 0;

  } // write block

  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename CONN_TYPE>
  int write_element_block(
      size_t blk_id,
      const std::string & name,
      const std::string & entity_desc,
      size_t num_elems,
      CONN_TYPE && element_conn)
  {
    return write_block<U>(
        blk_id, name, EX_ELEM_BLOCK, entity_desc, num_elems,
        std::forward<CONN_TYPE>(element_conn));
  }
  
  //============================================================================
  //! \brief write the global mapping a file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] ids  The ids to write
  //============================================================================
  template<typename U, typename V>
  int write_element_map(const V & ids) {
  
    using ex_index_t = U;
    std::vector<ex_index_t> mapping(ids.size());
    std::transform( ids.begin(), ids.end(), mapping.begin(),
        [](auto v){ return v+1; } );
    auto status = ex_put_elem_num_map(exoid_, mapping.data());
    if (status) return status;

    return 0;
  }

  //============================================================================
  //! \brief write the global mapping a file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] ids  The ids to write
  //============================================================================
  template<typename U, typename V>
  int write_node_map(const V & ids) {
  
    using ex_index_t = U;
    std::vector<ex_index_t> mapping(ids.size());
    std::transform( ids.begin(), ids.end(), mapping.begin(),
        [](auto v){ return v+1; } );
    auto status = ex_put_node_num_map(exoid_, mapping.data());
    if (status) return status;

    return 0;
  }
  
  //============================================================================
  //! \brief write the coordinates of the mesh from a file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] vertex_coord  the vertex coordinates
  //============================================================================
  template<typename V>
  int write_point_coords(const V & vertex_coord) {

    if (vertex_coord.empty())
      return 0;

    auto num_nodes = vertex_coord.size() / num_dims_;

    // exodus is kind enough to fetch the data in the real type we ask for
    auto status = ex_put_coord(
        exoid_,
        vertex_coord.data(),
        vertex_coord.data() + num_nodes,
        vertex_coord.data() + 2 * num_nodes);
    if (status) return status;
  
    return 0;
  }
  
  //============================================================================
  //! Write side sets
  //============================================================================
  template<typename U, typename V>
  int write_side_set(
      int_t ss_id,
      size_t num_sides,
      V && side_info_func)
  {
    // some type aliases
    using ex_index_t = U;

    // final storage
    auto status = ex_put_side_set_param(exoid_, ss_id, num_sides, 0);
    if (status) return status;
    
    if ( num_sides == 0 ) return 0;

    // pick out the sides
    std::vector<ex_index_t> elem_list;
    std::vector<ex_index_t> side_list;
    elem_list.reserve(num_sides);
    side_list.reserve(num_sides);

    ex_index_t elem, side;
    for (size_t s=0; s<num_sides; ++s) {
      std::forward<V>(side_info_func)(s, elem, side);
      elem_list.emplace_back(elem + 1);
      side_list.emplace_back(side + 1);
    }

    status = ex_put_side_set (exoid_, ss_id, elem_list.data(), side_list.data());
    if (status) return status;

    return 0;
  }

 
};

} // namepsace

#endif
