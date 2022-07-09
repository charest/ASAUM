// user includes
#include "csv_writer.hpp"
#include "comm/mpi_comm.hpp"
#include "lua/lua_utils.hpp"

namespace prl {

//==============================================================================
//! \brief Constructor with input
//==============================================================================
csv_writer_t::csv_writer_t(lua_result_t input, mpi_comm_t & comm) 
{
  auto is_root = comm.is_root();

  // significant figures
  if (!input["sigfigs"].empty()) {

    auto res = as_scalar<int_t>(input, "sigfigs", is_root);
    if (!res.exists) comm.exit(consts::failure);
    sigfigs_ = res.value;

    if (sigfigs_ < 0) {
      if (is_root) {
        std::cout << "'sigfigs' must be greater than or equal to 0, you"
          << " provided " << sigfigs_ << std::endl;
      }
      comm.exit(consts::failure);
    }
  
    width_ = sigfigs_ + 7; 
  }

  // seperator
  if (!input["seperator"].empty()) {
    auto res = as_scalar<std::string>(input, "seperator", is_root);
    if (!res.exists) comm.exit(consts::failure);
    sep_ = res.value;
  }

  // significant figures
  if (!input["width"].empty()) {

    auto res = as_scalar<int_t>(input, "width", is_root);
    if (!res.exists) comm.exit(consts::failure);
    width_ = res.value;

    if (width_ < 0) {
      if (is_root) {
        std::cout << "'width' must be greater than or equal to 0, you"
          << " provided " << width_ << std::endl;
      }
      comm.exit(consts::failure);
    }
  }

}

//==============================================================================
//! \brief Open a file
//==============================================================================
int csv_writer_t::open(const char* filename)
{
  file_.open(filename);
  file_.setf(std::ios::scientific);
  file_.precision(sigfigs_);
  col_ = 0;
  return !file_.good();
}
  
//==============================================================================
//! \brief Close a file
//==============================================================================
int csv_writer_t::close( void ) 
{   
  file_.close();
  return !file_.good();
}

//==============================================================================
//! \brief add a title
//==============================================================================
int csv_writer_t::comment(const char* str)
{
  // start a new line if we are not at the first col
  if (col_ > 0) file_ << std::endl;
  file_ << "# " << str << std::endl;
  col_ = 0;
  return !file_.good();
}

//==============================================================================
//! \brief add a new line
//==============================================================================
int csv_writer_t::new_row()
{
  file_ << std::endl;
  col_ = 0;
  return !file_.good();
}

    
} // namespace
