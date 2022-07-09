#ifndef CSV_WRITER_HPP
#define CSV_WRITER_HPP
 
#include "config.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>

namespace prl {

class lua_result_t;
class mpi_comm_t;

////////////////////////////////////////////////////////////////////////////////
/// csv writer class
////////////////////////////////////////////////////////////////////////////////
class csv_writer_t {
public:

  /// constructors
  csv_writer_t() = default;
  csv_writer_t(lua_result_t input, mpi_comm_t & comm);

  /// open and close
  int open(const char * filename);
  int close();

  int comment(const char * str);
  int new_row();
  
  //============================================================================
  /// Write a data array
  //============================================================================
  template< typename T>
  int write(const T & data)
  {
    if (col_) file_ << sep_;
    file_ << std::setw(width_) << data;
    col_++;
    return !file_.good();
  }  
  
  template< typename T>
  csv_writer_t & operator<<(const T & data)
  {
    write(data);
    return *this;
  }
  

private:

  //! \brief file pointer
  std::ofstream file_;
  
  //! significant figures
  int sigfigs_ = consts::digits;
  int width_ = consts::digits + 7; 

  //! seperator
  std::string sep_ = ",";

  //! column counter
  int col_ = 0;

};

} // namespace

#endif
