#ifndef EXO_READER_HPP
#define EXO_READER_HPP

#include "exo_types.hpp"
#include "shape.hpp"

#include <exodusII.h>

#include <algorithm>
#include <cstring>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace prl {


////////////////////////////////////////////////////////////////////////////////
/// an exodus reader class
////////////////////////////////////////////////////////////////////////////////
class exo_reader_t {
  int exoid_ = -1;
  int num_dim_ = 0;

public:

  int open(const std::string & name);
  int close();
  int read_params(ex_init_params &);
  bool is_int64();
  
  //============================================================================
  //! \brief read the block ids from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U, typename V>
  int read_block_ids(
      ex_entity_type obj_type,
      size_t num_blocks,
      std::vector<V> & ids)
  {
    // some type aliases
    using ex_index_t = U;
  
    // final storage
    ids.clear();
    ids.resize(num_blocks);
  
    if (num_blocks > 0) {
  
      // get the ids first
      std::vector<ex_index_t> block_ids(num_blocks);
      auto status = ex_get_ids(exoid_, obj_type, block_ids.data());
      if (status) return status;
      // now convert them
      std::copy(block_ids.begin(), block_ids.end(), ids.begin());
    }
  
    // now return them
    return 0;
  }

  
  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  int read_block_stats(
      ex_entity_id blk_id,
      ex_entity_type entity_type,
      exo_block_stats_t<U> & block_stats)
  {
    auto status = ex_get_block(
        exoid_,
        entity_type,
        blk_id,
        block_stats.elem_type,
        &block_stats.num_elem_this_blk,
        &block_stats.num_nodes_per_elem,
        &block_stats.num_edges_per_elem,
        &block_stats.num_faces_per_elem,
        &block_stats.num_attr);
    return status;
  }
  
  template<typename U>
  int read_block_stats(
      ex_entity_type entity_type,
      const std::vector<int_t> & elem_blk_ids,
      std::vector<exo_block_stats_t<U>> & block_stats)
  {
    auto n = elem_blk_ids.size();
    block_stats.clear();
    block_stats.resize(n);
    for (size_t i=0; i<n; ++i) {
      auto status = read_block_stats(
          elem_blk_ids[i],
          entity_type,
          block_stats[i]);
      if (status) return status;
    }
    return 0;
  }
  
  //============================================================================
  //! \brief read the element blocks from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  int read_block(
      ex_entity_id blk_id,
      ex_entity_type entity_type,
      const exo_block_stats_t<U> & stats,
      std::vector<int> & counts,
      std::vector<U> & indices,
      shape_t & block_info,
      bool &finished,
      U start = 0,
      U count = std::numeric_limits<U>::max())
  {
    finished = false;
    
    counts.clear();
    indices.clear();

    //------------------------------------------------------------------------
    // polygon data
    if (strcasecmp("nsided", stats.elem_type) == 0) {

      // the number of nodes per element is really the number of nodes
      // in the whole block
      auto num_nodes_this_blk = stats.num_nodes_per_elem;

      // get the number of nodes per element
      counts.resize(stats.num_elem_this_blk);
      auto status = ex_get_entity_count_per_polyhedra(
          exoid_, entity_type, blk_id, counts.data());
      if (status) return status;
      
      // read element definitions
      indices.resize(num_nodes_this_blk);
      status = ex_get_conn(
          exoid_, entity_type, blk_id, indices.data(), nullptr, nullptr);
      if (status) return status;

      block_info = shape_t::polygon;
      return 0;
    }
    //------------------------------------------------------------------------
    // polygon data
    else if (strcasecmp("nfaced", stats.elem_type) == 0) {

      // the number of faces per element is really the number of
      // faces in the whole block ( includes duplicate / overlapping
      // nodes )
      auto num_face_this_blk = stats.num_faces_per_elem;

      // get the number of nodes per element
      counts.resize(stats.num_elem_this_blk);
      auto status = ex_get_entity_count_per_polyhedra(
          exoid_, entity_type, blk_id, counts.data());
      if (status) return status;

      // read element definitions
      indices.resize(num_face_this_blk);
      status = ex_get_conn(
          exoid_, entity_type, blk_id, nullptr, nullptr, indices.data());
      if (status) return status;

      block_info = shape_t::polyhedron;
      return 0;

    }
    //------------------------------------------------------------------------
    // fixed element size
    else {

      auto end = std::min(start+count, stats.num_elem_this_blk);
      auto num_elem_this_blk = end - start;
      finished = (end == stats.num_elem_this_blk);

      // set the counts
      counts.resize(num_elem_this_blk);
      std::fill( counts.begin(), counts.end(), stats.num_nodes_per_elem );

      // read element definitions
      indices.resize(num_elem_this_blk * stats.num_nodes_per_elem);
      auto status = ex_get_partial_conn(
          exoid_,
          entity_type,
          blk_id,
          start+1,
          end-start,
          indices.data(),
          0,
          0);
      if (status) return status;

      // return element type
      if (
          strcasecmp("bar", stats.elem_type) == 0 ||
          strcasecmp("bar3", stats.elem_type) == 0)
      {
        block_info = shape_t::line;
      }
      else if (
          strcasecmp("tri", stats.elem_type) == 0 ||
          strcasecmp("tri3", stats.elem_type) == 0)
      {
        block_info = shape_t::tri;
      }
      else if (
          strcasecmp("quad", stats.elem_type) == 0 ||
          strcasecmp("quad4", stats.elem_type) == 0 ||
          strcasecmp("shell", stats.elem_type) == 0 ||
          strcasecmp("shell4", stats.elem_type) == 0)
      {
        block_info = shape_t::quad;
      }
      else if (
          strcasecmp("tet", stats.elem_type) == 0 ||
          strcasecmp("tetra", stats.elem_type) == 0 ||
          strcasecmp("tet4", stats.elem_type) == 0)
      {
        block_info = shape_t::tet;
      }
      else if (
          strcasecmp("hex", stats.elem_type) == 0 ||
          strcasecmp("hex8", stats.elem_type) == 0)
      {
        block_info = shape_t::hex;
      }
      else {
        std::cout << "Unknown exodus block type: " << stats.elem_type << std::endl;
        block_info = shape_t::unknown;
        return -1;
      }

    } // element type
    //------------------------------------------------------------------------
    
    return 0;
  }
 
  //============================================================================
  //! \brief read the global mapping from file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] ids  The ids to write
  //============================================================================
  template<
    typename U,
    typename V,
    bool Enabled = std::is_same<U, V>::value,
    typename = typename std::enable_if_t<Enabled>
  >
  int read_element_map(size_t num_elem, std::vector<V> & ids) {
    
    ids.clear();
    ids.resize(num_elem);
    
    auto status = ex_get_elem_num_map(exoid_, ids.data());
    if (status) return status;
    
    for ( auto & v : ids ) v--;
    
    return 0;
  }

  template<
    typename U,
    typename V,
    bool Enabled = !std::is_same<U, V>::value,
    typename std::enable_if_t<Enabled>* = nullptr
  >
  int read_element_map(size_t num_elem, std::vector<V> & ids) {
    using ex_index_t = U;
    std::vector<ex_index_t> mapping(num_elem);
    
    auto status = ex_get_elem_num_map(exoid_, mapping.data());
    if (status) return status;
    
    ids.clear();
    ids.reserve(num_elem);
    for ( auto v : mapping ) ids.emplace_back(v-1);
    
    return 0;
  }

  
  
  //============================================================================
  //! \brief read the global mapping from file.
  //! \param [in] exo_id  The exodus file id.
  //! \param [in] ids  The ids to write
  //============================================================================
  template<
    typename U,
    typename V,
    bool Enabled = std::is_same<U, V>::value,
    typename = typename std::enable_if_t<Enabled>
  >
  int read_node_map(size_t num_nodes, std::vector<V> & ids) {

    ids.clear();
    ids.resize(num_nodes);

    auto status = ex_get_node_num_map(exoid_, ids.data());
    if (status) return status;

    for ( auto & v : ids ) v--;
    
    return 0;
  }
  
  template<
    typename U,
    typename V,
    bool Enabled = !std::is_same<U, V>::value,
    typename std::enable_if_t<Enabled>* = nullptr
  >
  int read_node_map(size_t num_nodes, std::vector<V> & ids) {
    using ex_index_t = U;
    std::vector<ex_index_t> mapping(num_nodes);
    auto status = ex_get_node_num_map(exoid_, mapping.data());
    if (status) return status;

    ids.clear();
    ids.reserve(num_nodes);
    for ( auto v : mapping ) ids.emplace_back(v-1);
    
    return 0;
  }
  
  
  //============================================================================
  //! \brief read the coordinates of the mesh from a file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the vertex coordinates
  //============================================================================
  int read_point_coords(
      size_t num_dims,
      size_t num_nodes,
      std::vector<real_t> & vertex_coord);
  
  //============================================================================
  //! \brief read the coordinates of the mesh from a file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the vertex coordinates
  //============================================================================
  int read_point_coords(
      size_t num_dims,
      size_t start,
      size_t end,
      std::vector<real_t> & vertex_coord);


  //============================================================================
  //! \brief read the block ids from an exodus file.
  //! \param [in] exo_id  The exodus file id.
  //! \return the status of the file
  //============================================================================
  template<typename U>
  int read_block_params(
      ex_entity_type obj_type,
      const std::vector<U> & elem_blk_ids,
      std::vector<ex_elem_blk_parm> & elem_blk_parms)
  {
    elem_blk_parms.resize(elem_blk_ids.size());
    auto ndim = ex_inquire_int(exoid_, EX_INQ_DIM);

    ex_block block;
    size_t elem_ctr = 0;
    for (unsigned i=0; i<elem_blk_ids.size(); ++i) {

      block.id = elem_blk_ids[i];
      block.type = obj_type;
      auto status = ex_get_block_param (exoid_, &block);
      if (status) return status;
    
      elem_blk_parms[i].num_elem_in_blk = block.num_entry;
      elem_blk_parms[i].num_nodes_per_elem = block.num_nodes_per_entry;
      elem_blk_parms[i].num_attr = block.num_attribute;
      elem_blk_parms[i].elem_blk_id = block.id;    /* save id */
      elem_blk_parms[i].num_sides = 1;
   
      unsigned m=0;
      for (m=0; m < strlen(block.topology); m++)
        elem_blk_parms[i].elem_type[m] = std::toupper(block.topology[m]);
      elem_blk_parms[i].elem_type[m] = '\0';

    
      if (strncmp(elem_blk_parms[i].elem_type,"CIRCLE",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_CIRCLE;
	      /* set side set node stride */
        elem_blk_parms[i].num_nodes_per_side[0] = 1;
      }
      else if (strncmp(elem_blk_parms[i].elem_type,"SPHERE",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_SPHERE;
	      /* set side set node stride */
        elem_blk_parms[i].num_nodes_per_side[0] = 1;
      }
      else if (strncmp(elem_blk_parms[i].elem_type,"QUAD",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_QUAD;
	      /* determine side set node stride */
	      if (elem_blk_parms[i].num_nodes_per_elem == 4)
	        elem_blk_parms[i].num_nodes_per_side[0] = 2;
	      else if (elem_blk_parms[i].num_nodes_per_elem == 5)
	        elem_blk_parms[i].num_nodes_per_side[0] = 2;
	      else
	        elem_blk_parms[i].num_nodes_per_side[0] = 3;
      }
      else if (strncmp(elem_blk_parms[i].elem_type,"TRIANGLE",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_TRIANGLE;
	      /* set default side set node stride */
	      if (ndim == 2)  /* 2d TRIs */
	      {
	        if (elem_blk_parms[i].num_nodes_per_elem == 3)
	          elem_blk_parms[i].num_nodes_per_side[0] = 2;
	        else 
	          elem_blk_parms[i].num_nodes_per_side[0] = 3;
	      }
	      else if (ndim == 3)  /* 3d TRIs */
	      {
	        if (elem_blk_parms[i].num_nodes_per_elem == 3)
	          elem_blk_parms[i].num_nodes_per_side[0] = 3;
	        else 
	          elem_blk_parms[i].num_nodes_per_side[0] = 6;
	      }
      }
      else if (strncmp(elem_blk_parms[i].elem_type,"SHELL",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_SHELL;
	      /* determine side set node stride */
	      if (elem_blk_parms[i].num_nodes_per_elem == 2) /* KLUDGE for 2D Shells*/
	        elem_blk_parms[i].num_nodes_per_side[0] = 2;
	      else if (elem_blk_parms[i].num_nodes_per_elem == 4)
	        elem_blk_parms[i].num_nodes_per_side[0] = 4;
	      else
	        elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
      else if (strncmp(elem_blk_parms[i].elem_type,"HEX",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_HEX;
	      /* determine side set node stride */
	      if (elem_blk_parms[i].num_nodes_per_elem == 8)  /* 8-node bricks */
	        elem_blk_parms[i].num_nodes_per_side[0] = 4;
	      else if (elem_blk_parms[i].num_nodes_per_elem == 9)  /* 9-node bricks */
	        elem_blk_parms[i].num_nodes_per_side[0] = 4;
	      else if (elem_blk_parms[i].num_nodes_per_elem == 12)  /* HEXSHELLS */
	        elem_blk_parms[i].num_nodes_per_side[0] = 4;
	      else if (elem_blk_parms[i].num_nodes_per_elem == 27)  /* 27-node bricks */
	        elem_blk_parms[i].num_nodes_per_side[0] = 9;
	      else 
	        elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
      else if (strncmp(elem_blk_parms[i].elem_type,"TETRA",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_TETRA;
	      /* determine side set node stride */
	      if (elem_blk_parms[i].num_nodes_per_elem == 4)
	        elem_blk_parms[i].num_nodes_per_side[0] = 3;
	      else if (elem_blk_parms[i].num_nodes_per_elem == 8)
	        elem_blk_parms[i].num_nodes_per_side[0] = 4;
	      else
	        elem_blk_parms[i].num_nodes_per_side[0] = 6;
      }
      else if (strncmp(elem_blk_parms[i].elem_type,"WEDGE",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_WEDGE;
	      /* determine side set node stride */
	      if (elem_blk_parms[i].num_nodes_per_elem == 6)
	        elem_blk_parms[i].num_nodes_per_side[0] = 4;
	      else
	        elem_blk_parms[i].num_nodes_per_side[0] = 8;
            }
          else if (strncmp(elem_blk_parms[i].elem_type,"PYRAMID",3) == 0)
            {
	      elem_blk_parms[i].elem_type_val = EX_EL_PYRAMID;
	      /* determine side set node stride */
	      if (elem_blk_parms[i].num_nodes_per_elem == 5)
	        elem_blk_parms[i].num_nodes_per_side[0] = 4;
	      else
	        elem_blk_parms[i].num_nodes_per_side[0] = 8;
      }
      else if (strncmp(elem_blk_parms[i].elem_type,"BEAM",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_BEAM;
	      /* determine side set node stride */
	      if (elem_blk_parms[i].num_nodes_per_elem == 2)
	        elem_blk_parms[i].num_nodes_per_side[0] = 2;
	      else 
	        elem_blk_parms[i].num_nodes_per_side[0] = 3;
      }
      else if ( (strncmp(elem_blk_parms[i].elem_type,"TRUSS",3) == 0) ||
                (strncmp(elem_blk_parms[i].elem_type,"BAR",3) == 0) ||
                (strncmp(elem_blk_parms[i].elem_type,"EDGE",3) == 0) )
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_TRUSS;
	      /* determine side set node stride */
	      if (elem_blk_parms[i].num_nodes_per_elem == 2)
	        elem_blk_parms[i].num_nodes_per_side[0] = 2;
	      else 
	        elem_blk_parms[i].num_nodes_per_side[0] = 3;
      }
      else if (strncmp(elem_blk_parms[i].elem_type,"NULL",3) == 0)
      {
	      elem_blk_parms[i].elem_type_val = EX_EL_NULL_ELEMENT;
	      elem_blk_parms[i].num_nodes_per_side[0] = 0;
	      elem_blk_parms[i].num_elem_in_blk = 0;
      }
      else
      { /* unsupported element type; no problem if no sides specified for
	         this element block */
	      elem_blk_parms[i].elem_type_val = EX_EL_UNK;
	      elem_blk_parms[i].num_nodes_per_side[0] = 0;
      }

      elem_blk_parms[i].elem_blk_id = block.id;    /* save id */
      elem_ctr += elem_blk_parms[i].num_elem_in_blk;
      elem_blk_parms[i].elem_ctr = elem_ctr;      /* save elem number max */

    } // block

    return 0;

  }

  //============================================================================
  //! Read side set ids
  //============================================================================
  template<typename U, typename V>
  int read_side_set_ids(size_t num_side_sets, std::vector<V> & ids) {
    // some type aliases
    using ex_index_t = U;

    if (num_side_sets > 0) {

      // get the ids first
      std::vector<ex_index_t> ss_ids(num_side_sets);
      auto status = ex_get_side_set_ids(exoid_, ss_ids.data());
      if (status) return status;
      
      ids.clear();
      ids.resize(num_side_sets);

      // now convert them
      std::copy(ss_ids.begin(), ss_ids.end(), ids.begin());
    }

    // now return them
    return 0;
  }


  //============================================================================
  //! Read side set names
  //============================================================================
  int read_side_set_names(size_t num_side_sets, std::vector<std::string> & names);
  
  //============================================================================
  //! Read side set stats
  //============================================================================
  template<typename U>
  int read_side_set_stats(
      ex_entity_id ss_id,
      exo_side_set_stats_t<U> & side_set_stats
  ) {

    // get side set params
    auto status = ex_get_side_set_param(
        exoid_,
        ss_id,
        &side_set_stats.num_side_in_set,
        &side_set_stats.num_dist_fact_in_set);
    if (status) return status;

    // get side set connectivitiy lenght
    status = ex_get_side_set_node_list_len(
        exoid_,
        ss_id,
        &side_set_stats.side_set_node_list_len);
    if (status) return status;

    status = ex_get_set_param(exoid_, EX_SIDE_SET, ss_id, &side_set_stats.tot_num_ss_elem, 0);
    if (status) return status;

    return 0;
  }

  //============================================================================
  //! Read side sets
  //============================================================================
  template<typename U>
  int read_side_set(
      ex_entity_id ss_id,
      const exo_side_set_stats_t<U> & stats,
      std::vector<U> & side_set_elem_list,
      std::vector<U> & side_set_side_list,
      bool & finished,
      U start = 0,
      U count = std::numeric_limits<U>::max()
  ) {

    auto end = std::min(start+count, stats.tot_num_ss_elem);
    auto num_ss_elem = end - start;
    finished = (end == stats.tot_num_ss_elem);
    
    // get the actual connectivity
    side_set_elem_list.clear();
    side_set_side_list.clear();
    side_set_elem_list.resize(num_ss_elem);
    side_set_side_list.resize(num_ss_elem);
    auto status = ex_get_partial_side_set(
        exoid_,
        ss_id,
        start+1,
        end-start,
        side_set_elem_list.data(),
        side_set_side_list.data());
    if (status) return status;

    return 0;

  }
 

};

} // namepsace

#endif
