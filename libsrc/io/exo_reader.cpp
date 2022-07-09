#include "config.hpp"

#include "exo_reader.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Open and close a file
////////////////////////////////////////////////////////////////////////////////
int exo_reader_t::open(const std::string & name)
{ 
  // useful for debug
  //ex_opts(EX_ABORT | EX_VERBOSE);
    
  // size of floating point variables used in app.
  int app_word_size = sizeof(real_t);

  // size of floating point stored in name.
  int exo_word_size = 0;
  // the version number
  float version;

  // open the file
  exoid_ = ex_open(name.c_str(), EX_READ, &app_word_size, &exo_word_size, &version);
  if (exoid_ < 0) return consts::failure;
  
  //This sets the file to read IDs as 64 bit.  If the file does not have 
  //64 bit IDs, it should have no effect. 
  ex_set_int64_status(exoid_, EX_ALL_INT64_API);

  ex_init_params params;
  auto ret = ex_get_init_ext(exoid_, &params);
  if (ret) return ret;
  num_dim_ = params.num_dim;

  return consts::success;
}

////////////////////////////////////////////////////////////////////////////////
/// close a file
////////////////////////////////////////////////////////////////////////////////
int exo_reader_t::close()
{ return ex_close(exoid_); }
  
////////////////////////////////////////////////////////////////////////////////
/// read the mesh params
////////////////////////////////////////////////////////////////////////////////
int exo_reader_t::read_params(ex_init_params & params)
{
  // get the initialization parameters
  auto ret = ex_get_init_ext(exoid_, &params);
  return ret;
}

////////////////////////////////////////////////////////////////////////////////
/// check if this is a large integer support
////////////////////////////////////////////////////////////////////////////////
bool exo_reader_t::is_int64()
{ return (ex_int64_status(exoid_) & EX_IDS_INT64_API); }
  
//============================================================================
//! \brief read the coordinates of the mesh from a file.
//! \param [in] exo_id  The exodus file id.
//! \return the vertex coordinates
//============================================================================
int exo_reader_t::read_point_coords(
    size_t num_dims,
    size_t num_nodes,
    std::vector<real_t> & vertex_coord)
{

  // get the number of nodes
  if (num_nodes <= 0) return -1;
  
  // read nodes
  vertex_coord.clear();
  vertex_coord.resize(num_dims * num_nodes);

  // exodus is kind enough to fetch the data in the real type we ask for
  auto status = ex_get_coord(
      exoid_, vertex_coord.data(), vertex_coord.data() + num_nodes,
      vertex_coord.data() + 2 * num_nodes);
  if (status) return status;

  return 0;
}

//============================================================================
//! \brief read the coordinates of the mesh from a file.
//! \param [in] exo_id  The exodus file id.
//! \return the vertex coordinates
//============================================================================
int exo_reader_t::read_point_coords(
    size_t num_dims,
    size_t start,
    size_t end,
    std::vector<real_t> & vertex_coord) 
{

  // get the number of nodes
  auto num_nodes = end - start;
  if (num_nodes <= 0) return -1;

  // read nodes
  vertex_coord.clear();
  vertex_coord.resize(num_dims * num_nodes);

  // exodus is kind enough to fetch the data in the real type we ask for
  auto status = ex_get_partial_coord(
      exoid_,
      start+1,
      end-start,
      vertex_coord.data(),
      vertex_coord.data() + num_nodes,
      vertex_coord.data() + 2 * num_nodes);
  if (status) return status;

  return 0;
}
  
//============================================================================
//! Read side set names
//============================================================================
int exo_reader_t::read_side_set_names(
    size_t num_side_sets,
    std::vector<std::string> & names)
{

  if (num_side_sets > 0) {

    // get the ids first
    auto ss_names = new char*[num_side_sets];
    for (unsigned i=0; i<num_side_sets; ++i)
      ss_names[i] = new char[MAX_STR_LENGTH+1];
    auto status = ex_get_names (exoid_, EX_SIDE_SET, ss_names);
    if (status) return status;

    names.clear();
    names.resize(num_side_sets);

    for (unsigned i = 0; i < num_side_sets; i++){
      // if no label, use the id
      if ( strlen(ss_names[i]) != 0 )
        names [i] = ss_names[i];
    }
    
    // can clean up ss names
    for ( unsigned i=0; i<num_side_sets; ++i ) delete[] ss_names[i];
    delete[] ss_names;


  }

  // now return them
  return 0;
}
  

  
} // namespace
