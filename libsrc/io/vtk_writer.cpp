// user includes
#include "vtk_writer.hpp"
#include "comm/mpi_comm.hpp"
#include "lua/lua_utils.hpp"

namespace prl {

//==============================================================================
//! \brief Constructor with input
//==============================================================================
vtk_writer_t::vtk_writer_t(lua_result_t input, mpi_comm_t & comm) 
{
  auto is_root = comm.is_root();

  // Encoding type
  if (!input["encoding"].empty()) {
    
    auto res = validate<std::string>(
        input,
        "encoding",
        {"ascii", "binary", "zlib"},
      is_root);
    if (!res) comm.exit(consts::failure);
    auto str = input["encoding"].as<std::string>();

    if (str == "ascii") {
      encoding_ = encoding_t::ascii;
    }
    else if (str == "binary") {
      encoding_ = encoding_t::base64;
    }
    else if (str == "zlib") {
#ifndef HAVE_ZLIB
      if (is_root)
        std::cout << "Asking for 'zlib' compressed output, but "
          << "the code was built without zlib support.  Rebuild "
          << "with zlib to enable this option." << std::endl;
      comm.exit(consts::failure);
#endif
      encoding_ = encoding_t::zlib;
    }
  }

  // significant figures
  if (!input["sigfigs"].empty()) {

    auto res = as_scalar<int_t>(input, "sigfigs", is_root);
    if (!res.exists) comm.exit(consts::failure);
    sigfigs_ = res.value;

    if (sigfigs_ < 0) {
      if (is_root) {
        std::cout << "'sigfigs' must be greater than or equal to 0, you provided "
          << sigfigs_ << std::endl;
      }
      comm.exit(consts::failure);
    }
  }

}

//==============================================================================
//! \brief Open a file
//==============================================================================
int vtk_writer_t::open(const char* filename)
{
  file_.open(filename);
  file_ << "<?xml version=\"1.0\"?>" << std::endl;

  if (encoding_ == encoding_t::ascii) {
    file_.setf(std::ios::scientific);
    file_.precision(sigfigs_);
  }
  
  return !file_.good();
}
  
//==============================================================================
//! \brief Close a file
//==============================================================================
int vtk_writer_t::close( void ) 
{   
  file_ << "</VTKFile>";
  file_.close();
  return !file_.good();
}
    
//==============================================================================
//! \brief Start vtk header
//==============================================================================
void vtk_writer_t::start_vtkfile(const char * name)
{
  file_ << "<VTKFile type=\"" << name << "\" version=\"1.0\" byte_order=\"";

  if (isBigEndian()) file_ << "BigEndian";
  else               file_ << "LittleEndian";
  file_ << "\"" << " header_type=\"" << type_name<uint_t>() << "\"";
  if (encoding_ == encoding_t::zlib)
    file_ << " compressor=\"vtkZLibDataCompressor\"";
  file_ << ">" << std::endl;
}
  
//==============================================================================
//! \brief Start an unstructured mesh
//==============================================================================
int vtk_writer_t::start_unstructured() 
{   
  start_vtkfile("UnstructuredGrid");
  file_ << "<UnstructuredGrid>" << std::endl;
  return !file_.good();
}

//==============================================================================
//! \brief Stop an unstructured mesh
//==============================================================================
int vtk_writer_t::end_unstructured() 
{   
  file_ << "</UnstructuredGrid>" << std::endl;
  return !file_.good();
}

//==============================================================================
//! \brief Stop a structured mesh
//==============================================================================
int vtk_writer_t::end_structured() 
{   
  file_ << "</StructuredGrid>" << std::endl;
  return !file_.good();
}


//==============================================================================
//! \brief start an unsturctured piece
//==============================================================================
int vtk_writer_t::start_piece(std::size_t num_cells, std::size_t num_points) 
{   

  file_ << "<Piece NumberOfPoints=\"" << num_points << "\" "
    << "NumberOfCells=\"" << num_cells << "\">" << std::endl;

  return !file_.good();
}

//==============================================================================
/// End a piece
//==============================================================================
int vtk_writer_t::end_piece() {
  file_ << "</Piece>" << std::endl;
  return !file_.good();

}

//==============================================================================
/// start cell data
//==============================================================================
int vtk_writer_t::start_cell_data() 
{      
  file_ << "<CellData>" << std::endl;
  return !file_.good();
}

//==============================================================================
/// End cell data
//==============================================================================
int vtk_writer_t::end_cell_data() 
{      
  file_ << "</CellData>" << std::endl;
  return !file_.good();
}


//==============================================================================
/// Start point data
//==============================================================================
int vtk_writer_t::start_point_data() 
{   
  file_ << "<PointData>" << std::endl;
  return !file_.good();
}

//==============================================================================
/// End point data
//==============================================================================
int vtk_writer_t::end_point_data() 
{   
  file_ << "</PointData>" << std::endl;
  return !file_.good();
}
  
//==============================================================================
/// Start field data
//==============================================================================
int vtk_writer_t::start_field_data()
{   
  file_ << "<FieldData>" << std::endl;
  return !file_.good();
}

//==============================================================================
/// End field data
//==============================================================================
int vtk_writer_t::end_field_data()
{   
  file_ << "</FieldData>" << std::endl;
  return !file_.good();
}


} // namespace prl
