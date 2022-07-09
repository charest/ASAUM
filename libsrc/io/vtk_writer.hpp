#ifndef VTK_HPP
#define VTK_HPP

// user includes
#include "write_binary.hpp"
#include "write_zlib.hpp"
#include "utils/errors.hpp"
#include "utils/string_utils.hpp"

// system includes
#include <sstream>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace prl {

class lua_result_t;
class mpi_comm_t;
  
/// The element map.
enum class vtk_cell_type_t : std::uint8_t
{
  line = 3,
  triangle = 5,
  polygon = 7,
  quad = 9,
  tetra = 10,
  hexahedron = 12,
  wedge = 13,
  pyramid = 14,
  polyhedron = 42
};


////////////////////////////////////////////////////////////////////////////////
//! \brief A vtk writer class for legacy files.
////////////////////////////////////////////////////////////////////////////////
class vtk_writer_t {

public :

  /// Unsigned int type for vtk file
  using uint_t = std::uint64_t;

  /// the encoding type
  enum class encoding_t {
    ascii,
    binary, // broken cause needs to be appended to end of file
    base64,
    zlib
  };
  
  /// a header for writing compressed data
  struct header_t {
    uint_t blocks;
    uint_t blocksize;
    uint_t last_blocksize;
    uint_t compressed_blocksize;
  };

  //============================================================================
  //! \brief Determine the type string (int)
  //============================================================================
  template<
    typename T,
    typename = typename std::enable_if_t<std::is_integral<T>::value>
  >
  static const char * type_name() {
    if (std::is_signed<T>::value) {
      switch (sizeof(T)) {
        case 1: return "Int8";
        case 2: return "Int16";
        case 4: return "Int32";
        case 8: return "Int64";
      }
      THROW_ERROR("Unknown signed integer with size " << sizeof(T));
    }
    else {
      switch (sizeof(T)) {
        case 1: return "UInt8";
        case 2: return "UInt16";
        case 4: return "UInt32";
        case 8: return "UInt64";
      }
      THROW_ERROR("Unknown unsigned integer with size " << sizeof(T));
    }
    return "Unknown"; 
  }
  
  //============================================================================
  //! \brief Determine the type string (double)
  //============================================================================
  template<
    typename T,
    typename std::enable_if_t<std::is_floating_point<T>::value>* = nullptr
  >
  static const char * type_name() {
    switch (sizeof(T)) {
      case 4: return "Float32";
      case 8: return "Float64";
      THROW_ERROR("Unknown real with size " << sizeof(T));
    }
    return "Unknown"; 
  }
  
  //============================================================================
  //! \brief get the encoding string name
  //============================================================================
  static const char * encoding_format(encoding_t enc) {
    switch (enc) {
      case encoding_t::ascii:  return "ascii";
      case encoding_t::binary: return "binary";
      case encoding_t::base64: return "binary";
      case encoding_t::zlib:   return "binary";
    }
    THROW_ERROR("Unknown encoding type: " << static_cast<int>(enc));
  }

  //============================================================================
  //! \brief get the encoding string name
  //============================================================================
  static const char * encoding_type(encoding_t enc) {
    switch (enc) {
      case encoding_t::ascii:  return "";
      case encoding_t::binary: return "raw";
      case encoding_t::base64: return "base64";
      case encoding_t::zlib:   return "base64";
    }
    THROW_ERROR("Unknown encoding type: " << static_cast<int>(enc));
  }

  /// constructors
  vtk_writer_t() = default;
  vtk_writer_t(lua_result_t input, mpi_comm_t & comm);

  /// Open a file
  int open(const char* filename);

  /// Close a file
  int close( void );

private:

  void start_vtkfile(const char*name);

public:

  /// Start an unstructured mesh
  int start_unstructured();

  /// stop an unstructured mesh
  int end_unstructured();

  //============================================================================
  // \brief Start an structured mesh
  //============================================================================
  template<typename T>
  int start_structured(const std::vector<T> & dims)
  {
    auto num_dims = dims.size();

    start_vtkfile("StructuredGrid");

    file_ << "<StructuredGrid WholeExtent=\"";
    file_ << "0 " << dims[0];
    for (std::size_t d=1; d<num_dims; ++d) file_ << " 0 " << dims[d];
    for (std::size_t d=num_dims; d<3; ++d) file_ << " 0 0";

    file_ << "\">" << std::endl;
  
    return !file_.good();
  }

  // \brief Stop an structured mesh
  int end_structured();

  /// Start an unstructured piece
  int start_piece(std::size_t num_cells, std::size_t num_points);

  //============================================================================
  /// Start a structured piece
  //============================================================================
  template< typename U >
  int start_piece(const std::vector<U> & dims) {

    file_ << "<Piece Extent=\"";
    file_ << "0 " << dims[0];
    
    U nd = dims.size();
    for (U d=1; d<nd; ++d) file_ << " 0 " << dims[d];
    for (std::size_t d=nd; d<3; ++d) file_ << " 0 0";
    
    file_ << "\">" << std::endl;

    return !file_.good();

  }
  
  /// stop a piece
  int end_piece();
  
  /// start cell data
  int start_cell_data(); 
  /// stop cell data
  int end_cell_data(); 
  
  /// start point data
  int start_point_data(); 
  /// stop point data
  int end_point_data(); 
  
  /// start field data
  int start_field_data(); 
  /// stop field data
  int end_field_data(); 

  //============================================================================
  /// Write a data array
  //============================================================================
  template<typename C>
  int write_data_array(
      const C & data,
      std::size_t ncomp,
      const std::string & name = "",
      bool is_dimensional = false)
  {
    using value_type = typename C::value_type;

    auto nsize = data.size();
    auto nvals = nsize / ncomp;

    if (nsize != nvals * ncomp) {
      std::cout << "Field with name '" << name << "' has length of " 
        << nsize << " which is not evenly divisible by " << ncomp
        << std::endl;
      return -1;
    }

    // dirty hack for 1d/2d
    if (ncomp < 3 && is_dimensional) {
      std::vector<typename std::remove_cv<value_type>::type> newdat;
      newdat.reserve(nsize*3);
      for (size_t i=0; i<nvals; ++i) {
        for (size_t d=0; d<ncomp; ++d)
          newdat.emplace_back(data[i*ncomp + d]);
        for (int_t d=ncomp; d<3; ++d)
          newdat.emplace_back(value_type());
      }
      return write_data_array(newdat, 3, name);
    }

    std::string format = encoding_format(encoding_);

    file_ << "<DataArray type=\"" << type_name<value_type>() << "\"";
    if (!name.empty()) file_ << " Name=\"" << name << "\"";
    file_ << " NumberOfTuples=\"" << nvals <<  "\"";
    file_ << " NumberOfComponents=\"" << ncomp <<  "\"";
    file_ << " format=\"" << format << "\"";

    if (encoding_ != encoding_t::ascii) {
      auto enc = encoding_type(encoding_);
      file_ << " encoding=\"" << enc << "\"";
    }

    file_ << ">" << std::endl;

    // raw binary
    if (encoding_ == encoding_t::binary) {
      uint_t nbytes = nsize * sizeof(value_type);
      file_.write(
          reinterpret_cast<const char *>(&nbytes),
          sizeof(uint_t));
      file_.write(
          reinterpret_cast<const char *>(data.data()),
          nbytes);
    }
    // base64 encoded binary
    else if (encoding_ == encoding_t::base64) {
      uint_t nbytes = nsize * sizeof(value_type);
      auto buf = base64_encode(
          reinterpret_cast<const byte_t*>(&nbytes),
          sizeof(uint_t));
      file_ << buf;
      buf = base64_encode(
          reinterpret_cast<const byte_t*>(data.data()),
          nbytes);
      file_ << buf;
    }
    // compressed binary
    else if (encoding_ == encoding_t::zlib) {
      uint_t nbytes = nsize * sizeof(value_type);
      auto defdata = zlib_deflate(
          reinterpret_cast<const byte_t*>(data.data()),
          nbytes);
      header_t header;
      header.blocks = 1;
      header.blocksize = nbytes;
      header.last_blocksize = nbytes;
      header.compressed_blocksize = defdata.size();
      auto buf = base64_encode(
          reinterpret_cast<const byte_t*>(&header),
          sizeof(header_t));
      file_ << buf;
      buf = base64_encode(defdata.data(), defdata.size());
      file_ << buf;
    }
    // asci
    else {
      for (size_t i=0; i<nsize; ++i )
        file_ << data[i] << " ";
    }

    file_ << std::endl << "</DataArray>" << std::endl;

    return !file_.good();
               
  }  
  
  //============================================================================
  /// Write a data array
  //============================================================================
  template<typename T>
  int write_data_array(
      const T & data,
      const std::string & name = "")
  {
    std::string format = encoding_format(encoding_);

    file_ << "<DataArray type=\"" << type_name<T>() << "\"";
    if (!name.empty()) file_ << " Name=\"" << name << "\"";
    file_ << " NumberOfTuples=\"1\"";
    file_ << " format=\"" << format << "\"";

    if (encoding_ != encoding_t::ascii) {
      auto enc = encoding_type(encoding_);
      file_ << " encoding=\"" << enc << "\"";
    }
    
    file_ << ">" << std::endl;
    
    // raw binary
    if (encoding_ == encoding_t::binary) {
      uint_t nbytes = sizeof(T);
      file_.write(
          reinterpret_cast<const char *>(&nbytes),
          sizeof(uint_t));
      file_.write(
          reinterpret_cast<const char *>(&data),
          nbytes);
    }
    // base64 encoded binary
    else if (encoding_ == encoding_t::base64) {
      uint_t nbytes = sizeof(T);
      auto buf = base64_encode(
          reinterpret_cast<const byte_t*>(&nbytes),
          sizeof(uint_t));
      file_ << buf;
      buf = base64_encode(
          reinterpret_cast<const byte_t*>(&data),
          nbytes);
      file_ << buf;
    }
    // compressed binary
    else if (encoding_ == encoding_t::zlib) {
      uint_t nbytes = sizeof(T);
      auto defdata = zlib_deflate(
          reinterpret_cast<const byte_t*>(&data),
          nbytes);
      header_t header;
      header.blocks = 1;
      header.blocksize = nbytes;
      header.last_blocksize = nbytes;
      header.compressed_blocksize = defdata.size();
      auto buf = base64_encode(
          reinterpret_cast<const byte_t*>(&header),
          sizeof(header_t));
      file_ << buf;
      buf = base64_encode(defdata.data(), defdata.size());
      file_ << buf;
    }
    // ascii
    else {
      file_ << data;
    }

    file_ << std::endl <<  "</DataArray>" << std::endl;

    return !file_.good();
               
  }  


  //============================================================================
  /// Write point coordinates
  //============================================================================
  template<typename T>
  int write_points( const T & data, std::size_t num_dims )
  {

    file_ << "<Points>" << std::endl;
    write_data_array(data, num_dims, "Points", true);
    file_ << "</Points>" << std::endl;

    return !file_.good();
               
  }  



  /*! *****************************************************************
   * \brief Write connectivity information, i.e. cell to vertex
   *        connectivity.
   *
   * The number of vertices for each cell is ascertained from the 
   * cell type.
   *
   * \param [in] data  The connectivity data to write.  This is a 
   *                   flat array with the vertex ids listed for
   *                   each cell.
   * \param [in] cell_type The array of cell type flags for each cell.
   * \tparam C The container class the data is stored in.
   * \tparam T  The type of data stored in the container.
   * \tparam Args The rest of the args in the container.
   * \return 0 for success, 1 otherwise.   
   ********************************************************************/
  template< 
    typename OffsetsC,
    typename IndicesC,
    typename TypesC
  >
  int write_elements(
      const OffsetsC & offsets,
      const IndicesC & indices,
      const TypesC & types)
  {
  
    if (offsets.empty()) return !file_.good();

    auto offsets_ref = make_array_ref(offsets.data()+1, offsets.size()-1);

    file_ << "<Cells>" << std::endl;
    write_data_array(indices, 1, "connectivity");
    write_data_array(offsets_ref, 1, "offsets");
    write_data_array(types, 1, "types");
    file_ << "</Cells>" << std::endl;

    // check for write errors
    return !file_.good();
  
  }  


  template< 
    typename C2F_OffsetsC,
    typename C2F_IndicesC,
    typename C2V_OffsetsC,
    typename C2V_IndicesC,
    typename F2V_OffsetsC,
    typename F2V_IndicesC,
    typename TypesC
  >
  int write_elements(
      const C2F_OffsetsC & cell2face_offsets,
      const C2F_IndicesC & cell2face_indices,
      const C2V_OffsetsC & cell2vert_offsets,
      const C2V_IndicesC & cell2vert_indices,
      const F2V_OffsetsC & face2vert_offsets,
      const F2V_IndicesC & face2vert_indices,
      const TypesC & types)
  {
    
    if (cell2face_offsets.empty()) return !file_.good();

    file_ << "<Cells>" << std::endl;
    
    size_t ncell = cell2face_offsets.size() - 1;
    auto cell2vert_offsets_ref = make_array_ref( cell2vert_offsets.data()+1, ncell );

    write_data_array(cell2vert_indices, 1, "connectivity");
    write_data_array(cell2vert_offsets_ref, 1, "offsets");
    write_data_array(types, 1, "types");
    

    size_t cnt = 0;
    for (size_t c=0; c<ncell; ++c) {
      //--- cell -> face
      // length
      cnt += 1;
      //--- face -> vert
      for (auto f=cell2face_offsets[c]; f<cell2face_offsets[c+1]; ++f) {
        auto fid = cell2face_indices[f];
        // length + verts
        cnt += 1 + face2vert_offsets[fid+1] - face2vert_offsets[fid];
      }
    }
    
    using offset_t = typename F2V_OffsetsC::value_type;
    std::vector<offset_t> offsets;
    offsets.reserve(ncell);
    
    using index_t = typename F2V_IndicesC::value_type;
    std::vector<index_t> indices;
    indices.reserve(cnt); 

    for (size_t c=0; c<ncell; ++c) {
      //--- cell -> face
      auto fstart = cell2face_offsets[c];
      auto fend = cell2face_offsets[c+1];
      indices.emplace_back(fend - fstart);
      for (auto f=fstart; f<fend; ++f) {
        auto fid = cell2face_indices[f];
        //--- face -> vert
        auto vstart = face2vert_offsets[fid];
        auto vend = face2vert_offsets[fid+1];
        indices.emplace_back(vend - vstart);
        for (auto v=vstart; v<vend; ++v)
          indices.emplace_back( face2vert_indices[v] );
      }
      //--- offsets
      offsets.emplace_back(indices.size());
    }

    write_data_array(indices, 1, "faces");
    write_data_array(offsets, 1, "faceoffsets");
    
    file_ << "</Cells>" << std::endl;

    // check for write errors
    return !file_.good();
  
  }  


private:
  
  //============================================================================
  /// Write a mutliblock field
  //============================================================================
  void write_multiblock_field_data()
  {}

  template<typename T, typename...ARGS>
  void write_multiblock_field_data(
      const char * name,
      const T & data,
      const ARGS&... args)
  {
    write_data_array(data, name);
    return write_multiblock_field_data(args...);
  }

public:

  //============================================================================
  /// Write a mutliblock dataset
  //============================================================================
  template<typename T, typename...ARGS>
  int write_multiblock(
      const char * filename,
      const char * prefix,
      const std::vector<std::string> & postfixes,
      std::size_t ndigits,
      const std::vector<T> & block_offsets,
      const ARGS&... args)
  {

    std::size_t num_groups = block_offsets.empty() ?
      0 : block_offsets.size()-1;

    open(filename);

    start_vtkfile("vtkMultiBlockDataSet");
    file_ << "<vtkMultiBlockDataSet>" << std::endl;

    for (std::size_t g=0; g<num_groups; ++g) {
      for  (auto b=block_offsets[g]; b<block_offsets[g+1]; ++b) {
        std::stringstream ss;
        ss << prefix << zero_padded(b, ndigits) << postfixes[b];
        file_ << "<DataSet group=\"" << g  << "\" index=\"" << b << "\" file=\"";
        file_ << ss.str() << "\"/>" <<  std::endl;
      }
    }
    
    file_ << "</vtkMultiBlockDataSet>" << std::endl;


    file_ << "<FieldData>" << std::endl;
    write_multiblock_field_data(args...); 
    file_ << "</FieldData>" << std::endl;

    close();

    auto ret = file_.good();

    // check for errors
    return !ret;
  }


private :

  //! \brief file pointer
  std::ofstream file_;
  
  //! significant figures
  int_t sigfigs_ = consts::digits;

  //! \brief type of encoding
  encoding_t encoding_ = encoding_t::base64;

};


} // namespace

#endif
