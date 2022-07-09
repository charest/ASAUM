#include "comm_info.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// The dense array class
////////////////////////////////////////////////////////////////////////////////
size_t dense_info_t::pack(
    const array_ref<const int_t> & ids,
    byte_t * buf) const
{
  auto ptr = static_cast<byte_t*>(data);
  for (auto id : ids) {
    auto src = ptr + data_size*id;
    std::memcpy(buf, src, data_size);
    buf += data_size;
  } // i
  return ids.size() * data_size;
}

size_t dense_info_t::unpack(
    const array_ref<const int_t> & ids,
    const byte_t * buf)
{
  auto ptr = static_cast<byte_t*>(data);
  for (auto id : ids) {
    auto dst = ptr + data_size*id;
    std::memcpy(dst, buf, data_size);
    buf += data_size;
  } // i
  return ids.size() * data_size;
}
  
void dense_info_t::copy(
    const comm_info_t & src,
    const array_ref<const int_t> & src_ids,
    const array_ref<const int_t> & dst_ids)
{
  auto dst_ptr = static_cast<byte_t*>(data);
  auto src_ptr = static_cast<byte_t*>(src.data);
  
  if (data_size != src.data_size)
    THROW_ERROR("Source and destination communication queues have "
        << "differing field data sizes.");

  auto num = dst_ids.size();
  for (decltype(num) i=0; i<num; ++i) {
    auto dst_id = dst_ids[i];
    auto src_id = src_ids[i];
    auto dst = dst_ptr + data_size*dst_id;
    auto src = src_ptr + data_size*src_id;
    std::memcpy(dst, src, data_size);
  } // i
}

////////////////////////////////////////////////////////////////////////////////
/// The crs array class
////////////////////////////////////////////////////////////////////////////////
size_t crs_info_t::count(const array_ref<const int_t> & ids) const
{ 
  size_t sz = 0;
  if (counts) {
    for (auto i : ids)
      sz += counts[i];
  }
  else {
    for (auto i : ids)
      sz += offsets[i+1] - offsets[i];
  }
  return sz*data_size;
}

size_t crs_info_t::pack(
    const array_ref<const int_t> & ids,
    byte_t * buf) const
{
  auto ptr = static_cast<byte_t*>(data);
  auto start = buf;
  for (auto i : ids) {
    auto off = offsets[i];
    auto len = counts ? counts[i] : offsets[i+1] - off;
    auto nbytes = data_size * len;
    auto src = ptr + data_size*off;
    std::memcpy(buf, src, nbytes);
    buf += nbytes;
  } // i
  return std::distance(start, buf);
}

size_t crs_info_t::unpack(
    const array_ref<const int_t> & ids,
    const byte_t * buf)
{
  auto ptr = static_cast<byte_t*>(data);
  auto start = buf;
  for (auto i : ids) {
    auto off = offsets[i];
    auto len = counts ? counts[i] : offsets[i+1] - off;
    auto nbytes = data_size * len;
    auto dst = ptr + data_size*off;
    std::memcpy(dst, buf, nbytes);
    buf += nbytes;
  } // i
  return std::distance(start, buf);
}

void crs_info_t::copy(
    const comm_info_t & src,
    const array_ref<const int_t> & src_ids,
    const array_ref<const int_t> & dst_ids)
{
  auto dst_ptr = static_cast<byte_t*>(data);
  auto src_ptr = static_cast<byte_t*>(src.data);
  
  if (data_size != src.data_size)
    THROW_ERROR("Source and destination communication queues have "
        << "differing field data sizes.");

  auto src_crs = dynamic_cast<const crs_info_t*>(&src);
  if (!src_crs) THROW_ERROR("Source and destination are not the same types!");
  auto & src_offsets = src_crs->offsets;
  auto & src_counts = src_crs->counts;

  auto num = dst_ids.size();
  for (decltype(num) i=0; i<num; ++i) {
    auto dst_id = dst_ids[i];
    auto dst_off = offsets[dst_id];
    auto dst_len = counts ? counts[dst_id] : offsets[dst_id+1] - dst_off;
    
    auto src_id = src_ids[i];
    auto src_off = src_offsets[src_id];
    auto src_len = src_counts ? src_counts[src_id] : src_offsets[src_id+1] - src_off;

    assert(dst_len == src_len && "mismatch in crs_info_t::copy()");
    UNUSED(src_len);
    auto nbytes = data_size * dst_len;

    auto dst = dst_ptr + data_size*dst_off;
    auto src = src_ptr + data_size*src_off;
    std::memcpy(dst, src, nbytes);
  } // i
}

} // namespace
