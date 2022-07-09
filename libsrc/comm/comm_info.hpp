#ifndef COMM_INFO_HPP
#define COMM_INFO_HPP

#include "config.hpp"
#include "utils/array_ref.hpp"
#include "utils/errors.hpp"

#include <cstring>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Base class
////////////////////////////////////////////////////////////////////////////////
struct comm_info_t {
  int_t data_size = 0;
  void * data = nullptr;

  comm_info_t(int_t sz, void * d) : data_size(sz), data(d) {}
  
  virtual size_t count(const array_ref<const int_t> & ids) const = 0;
  virtual size_t pack(
      const array_ref<const int_t> & ids,
      byte_t * buf) const = 0;
  virtual size_t unpack(
      const array_ref<const int_t> & ids,
      const byte_t * buf) = 0;
  virtual void copy(
      const comm_info_t & src,
      const array_ref<const int_t> & src_ids,
      const array_ref<const int_t> & dst_ids) = 0;
  
  virtual ~comm_info_t() = default;
};

////////////////////////////////////////////////////////////////////////////////
/// The dense array class
////////////////////////////////////////////////////////////////////////////////
struct dense_info_t : public comm_info_t
{
  dense_info_t(int_t sz, void * d) : comm_info_t(sz, d) {}

  size_t count(const array_ref<const int_t> & ids) const override
  { return ids.size()*data_size; }

  size_t pack(
      const array_ref<const int_t> & ids,
      byte_t * buf) const override;
  
  size_t unpack(
      const array_ref<const int_t> & ids,
      const byte_t * buf) override;
    
  void copy(
      const comm_info_t & src,
      const array_ref<const int_t> & src_ids,
      const array_ref<const int_t> & dst_ids) override;

};

////////////////////////////////////////////////////////////////////////////////
/// The crs array class
////////////////////////////////////////////////////////////////////////////////
struct crs_info_t : public comm_info_t
{
  int_t * offsets = nullptr;
  int_t * counts = nullptr;

  crs_info_t(int_t sz, void * d, int_t * off) :
    comm_info_t(sz, d), offsets(off) {}
  
  crs_info_t(int_t sz, void * d, int_t * off, int_t * cnt) :
    comm_info_t(sz, d), offsets(off), counts(cnt) {}
  
  size_t count(const array_ref<const int_t> & ids) const override;
  
  size_t pack(
      const array_ref<const int_t> & ids,
      byte_t * buf) const override;
  
  size_t unpack(
      const array_ref<const int_t> & ids,
      const byte_t * buf) override;
  
  void copy(
      const comm_info_t & src,
      const array_ref<const int_t> & src_ids,
      const array_ref<const int_t> & dst_ids) override;
};

} // namespace

#endif
