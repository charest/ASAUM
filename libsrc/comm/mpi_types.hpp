#ifndef MPI_TYPES_HPP
#define MPI_TYPES_HPP

#include "comm_mutex.hpp"

#include <mpi.h>

#include <complex>
#include <list>
#include <memory>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Registry to keep track of data types
////////////////////////////////////////////////////////////////////////////////
struct mpi_type_registry_t {
  
  using value_type = std::unique_ptr<MPI_Datatype, void(*)(MPI_Datatype*)>;

  std::list<value_type> data_types;

  static mpi_type_registry_t& instance()
  {
    static mpi_type_registry_t s;
    return s;
  }
};

////////////////////////////////////////////////////////////////////////////////
/// MPI type traits helpers
////////////////////////////////////////////////////////////////////////////////
template<typename T, typename U = void>
struct mpi_type_t {

  static MPI_Datatype value() {

    // first create the type
    static MPI_Datatype * ret = [] {

      mpi_type_registry_t::value_type data_type(
          new MPI_Datatype(),
          [](MPI_Datatype * p) {
            comm_mutex_t::instance().lock();
            MPI_Type_free(p);
            comm_mutex_t::instance().unlock();
            delete p;
          }
      );

      comm_mutex_t::instance().lock();
      MPI_Type_contiguous(sizeof(T), MPI_BYTE, data_type.get());
      MPI_Type_commit(data_type.get());
      comm_mutex_t::instance().unlock();
      
      auto & reg = mpi_type_registry_t::instance().data_types;
      reg.emplace_back(std::move(data_type));
      
      return reg.back().get();
    }();
  
    // return thre resut
    return *ret;

  } // value

};

template<>
struct mpi_type_t<char>
{
  static MPI_Datatype value() { return MPI_CHAR; }
};

template<>
struct mpi_type_t<short>
{
  static MPI_Datatype value() { return MPI_SHORT; }
};

template<>
struct mpi_type_t<int>
{
  static MPI_Datatype value() { return MPI_INT; }
};

template<>
struct mpi_type_t<long>
{
  static MPI_Datatype value() { return MPI_LONG; }
};

template<>
struct mpi_type_t<long long>
{
  static MPI_Datatype value() { return MPI_LONG_LONG; }
};

template<>
struct mpi_type_t<signed char>
{
  static MPI_Datatype value() { return MPI_SIGNED_CHAR; }
};

template<>
struct mpi_type_t<unsigned char>
{
  static MPI_Datatype value() { return MPI_UNSIGNED_CHAR; }
};

template<>
struct mpi_type_t<unsigned short>
{
  static MPI_Datatype value() { return MPI_UNSIGNED_SHORT; }
};

template<>
struct mpi_type_t<unsigned int>
{
  static MPI_Datatype value() { return MPI_UNSIGNED; }
};

template<>
struct mpi_type_t<unsigned long>
{
  static MPI_Datatype value() { return MPI_UNSIGNED_LONG; }
};

template<>
struct mpi_type_t<unsigned long long>
{
  static MPI_Datatype value() { return MPI_UNSIGNED_LONG_LONG; }
};

template<>
struct mpi_type_t<float>
{
  static MPI_Datatype value() { return MPI_FLOAT; }
};

template<>
struct mpi_type_t<double>
{
  static MPI_Datatype value() { return MPI_DOUBLE; }
};

template<>
struct mpi_type_t<long double>
{
  static MPI_Datatype value() { return MPI_LONG_DOUBLE; }
};

template<>
struct mpi_type_t<wchar_t>
{
  static MPI_Datatype value() { return MPI_WCHAR; }
};

template<>
struct mpi_type_t<bool>
{
  static MPI_Datatype value() { return MPI_CXX_BOOL; }
};

template<>
struct mpi_type_t<std::complex<float>>
{
  static MPI_Datatype value() { return MPI_CXX_FLOAT_COMPLEX; }
};

template<>
struct mpi_type_t<std::complex<double>>
{
  static MPI_Datatype value() { return MPI_CXX_DOUBLE_COMPLEX; }
};

template<>
struct mpi_type_t<std::complex<long double>>
{
  static MPI_Datatype value() { return MPI_CXX_LONG_DOUBLE_COMPLEX; }
};


} // namepsace

#endif
