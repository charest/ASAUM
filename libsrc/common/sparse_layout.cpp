#include "sparse_layout.hpp"

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// Setup a sparse layout
////////////////////////////////////////////////////////////////////////////////
void sparse_layout_t::setup(
    int_t ncols,
    const std::vector<std::vector<int_t>> & ids,
    std::vector<int_t> & off,
    std::vector<int_t> & cnts,
    std::vector<int_t> & ind)
{
  num_columns = ncols;
  bandwidth = 0;
  row_offsets = nullptr;
  row_counts = nullptr;
  col_ids = nullptr;
    
  num_rows = ids.size();

  off.clear();
  cnts.clear();
  ind.clear();
    
  // max bandwidth
  for (const auto & row : ids)
    bandwidth = std::max<int_t>(row.size(), bandwidth);


  //----------------------------------
  // Itpack
  if (type == sparse_layout_type_t::itpack) {

    cnts.reserve(num_rows);
    ind.reserve(bandwidth*num_rows);

    for (const auto & row : ids) {
      auto n = row.size();
      for (size_t j=0; j<n; ++j)
        ind.emplace_back(row[j]);
      for (size_t j=n; j<bandwidth; ++j)
        ind.emplace_back(-1);
      cnts.emplace_back(n);
    }

    col_ids = ind.data();
    row_counts = cnts.data();

  }
  //----------------------------------
  // crs
  else if (type == sparse_layout_type_t::crs) {

    // count entries
    size_t len = 0;
    for (const auto & row : ids)
      len += row.size();

    // resize
    ind.reserve(len);
    off.reserve(num_rows + 1);
    cnts.reserve(num_rows);

    // fill
    off.emplace_back(0);
    for (const auto & row : ids) {
      auto n = row.size();
      for (size_t j=0; j<n; ++j)
        ind.emplace_back(row[j]);
      cnts.emplace_back(n);
      off.emplace_back(ind.size());
    }

    col_ids = ind.data();
    row_offsets = off.data();
    row_counts = cnts.data();

  }

} // make

////////////////////////////////////////////////////////////////////////////////
/// Setup a sparse layout
////////////////////////////////////////////////////////////////////////////////
void sparse_layout_t::setup(
    size_t tot_cols,
    std::vector<char> & active,
    std::vector<int_t> & off,
    std::vector<int_t> & cnts,
    std::vector<int_t> & ind)
{
  num_columns = tot_cols;
  bandwidth = 0;
  row_offsets = nullptr;
  row_counts = nullptr;
  col_ids = nullptr;
    
  num_rows = active.size() / tot_cols;

  off.clear();
  cnts.clear();
  ind.clear();
  
  // max bandwidth
  // count entries
  size_t len = 0;
  for (size_t i=0; i<num_rows; ++i) {
    size_t row_size{0};
    for (size_t j=0; j<num_columns; ++j)
      row_size += (active[i*tot_cols + j]);
    len += row_size;
    bandwidth = std::max<int_t>(row_size, bandwidth);
  }

    
  //----------------------------------
  // Itpack
  if (type == sparse_layout_type_t::itpack) {

    cnts.reserve(num_rows);
    ind.reserve(bandwidth*num_rows);

    for (size_t i=0; i<num_rows; ++i) {
      size_t row_cnt = 0;
      for (size_t j=0; j<num_columns; ++j) {
        if (active[i*tot_cols + j]) {
          ind.emplace_back(j);
          ++row_cnt;
        }
        else {
          ind.emplace_back(-1);
        }
        cnts.emplace_back(row_cnt);
      }
    }

    col_ids = ind.data();
    row_counts = cnts.data();

  }
  //----------------------------------
  // crs
  else if (type == sparse_layout_type_t::crs) {

    // resize
    ind.reserve(len);
    off.reserve(num_rows + 1);
    cnts.reserve(num_rows);

    // fill
    off.emplace_back(0);
    for (size_t i=0; i<num_rows; ++i) {
      auto start = ind.size();
      for (size_t j=0; j<num_columns; ++j)
        if (active[i*tot_cols + j])
          ind.emplace_back(j);
      auto end  = ind.size();
      off.emplace_back(end);
      cnts.emplace_back(end-start);
    }

    col_ids = ind.data();
    row_offsets = off.data();
    row_counts = cnts.data();

  }

} // make
  
  
//////////////////////////////////////////////////////////////////////////////
/// Re-Setup sparse representation
//////////////////////////////////////////////////////////////////////////////
void sparse_layout_t::reconfigure(
    size_t tot_cols,
    size_t max_bandwidth,
    const std::vector<int_t> & ncols,
    std::vector<int_t> & off,
    std::vector<int_t> & cnts,
    std::vector<int_t> & ind)
{
  //----------------------------------
  // Itpack
  if (type == sparse_layout_type_t::itpack) {
    resize<int_t>(tot_cols, max_bandwidth, ncols, ind, -1);

    for (size_t i=0; i<num_rows; ++i)
      cnts[i] = ncols[i];

    row_counts = cnts.data();
    col_ids = ind.data();
  }
  //----------------------------------
  // crs
  else if (type == sparse_layout_type_t::crs) {

    resize<int_t>(tot_cols, max_bandwidth, ncols, ind, -1);
    
    for (size_t i=0; i<num_rows; ++i) {
      auto n = ncols[i];
      cnts[i] = n;
      off[i+1] = off[i] + ncols[i];
    }
    
    row_offsets = off.data();
    row_counts = cnts.data();
    col_ids = ind.data();
  }
  
  // set some general stuff
  num_columns = tot_cols;
  bandwidth = max_bandwidth;


}

//////////////////////////////////////////////////////////////////////////////
/// Re-Setup sparse representation
//////////////////////////////////////////////////////////////////////////////
void sparse_layout_t::reset_columns(std::vector<char> & active)
{

  auto tot_cols = active.size() / num_rows;
  int_t num_cols = num_columns;

  //----------------------------------
  // Itpack
  if (type == sparse_layout_type_t::itpack) {

    for (size_t i=0; i<num_rows; ++i) {
      auto start = bandwidth*i;
      size_t cnt = 0;
      for (int_t j=0; j<num_cols; ++j) {
        if (active[i*tot_cols + j]) {
          assert(cnt<bandwidth && "itpack out of bounds");
          col_ids[start + cnt] = j;
          cnt++;
        }
      }
    }

  }
  //----------------------------------
  // crs
  else if (type == sparse_layout_type_t::crs) {
    
    for (size_t i=0; i<num_rows; ++i) {
      auto start = row_offsets[i];
      int_t cnt = 0;
      for (int_t j=0; j<num_cols; ++j) {
        if (active[i*tot_cols + j]) {
          assert(cnt<row_offsets[i+1]-row_offsets[i] && "crs out of bounds");
          col_ids[start + cnt] = j;
          cnt++;
        }
      }
    }

  }

}

//////////////////////////////////////////////////////////////////////////////
/// Check if a resize is needed
//////////////////////////////////////////////////////////////////////////////
bool sparse_layout_t::needs_resize(
    size_t tot_cols,
    size_t max_bandwidth,
    const std::vector<int_t> & ncols,
    bool shrink) const
{
  //--- dense storage
  if (type == sparse_layout_type_t::dense) {
    return shrink ? tot_cols!=num_columns : tot_cols>num_columns;
  }
  //--- max bandwidth
  else if (type == sparse_layout_type_t::itpack) {
      return shrink ? max_bandwidth!=bandwidth : max_bandwidth>bandwidth;
  }
  //--- compressed row
  else if (type == sparse_layout_type_t::crs) {
    if (shrink) {
      for (size_t i=0; i<num_rows; ++i)
        if (ncols[i]!=row_counts[i])
          return true;
    }
    else {
      for (size_t i=0; i<num_rows; ++i) {
        auto sz = row_offsets[i+1] - row_offsets[i];
        if (ncols[i] > sz)
          return true;
      }
    }
  }
  //--- no resize needed
  return false;
}

} // namesapce
