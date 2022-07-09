#ifndef SECTOR_HPP
#define SECTOR_HPP

#include "config.hpp"
#include "math/morton.hpp"

#include <ostream>
#include <unordered_map>
#include <vector>

namespace prl {

struct range_ijk_t;
struct sector_t;

extern const std::vector<sector_t> sector_permutations[3];

////////////////////////////////////////////////////////////////////////////////
/// Sector class
////////////////////////////////////////////////////////////////////////////////
struct sector_t {
  int_t sector_[3] = {0, 0, 0};
  
  sector_t() = default;
  sector_t(int_t i, int_t j, int_t k) : sector_{i, j, k} {}

  sector_t(bool is_hi, int_t dim) 
  { sector_[dim] = is_hi ? 1 : -1; }

  void set(bool is_hi, int_t dim) 
  { 
    sector_[0] = 0;
    sector_[1] = 0;
    sector_[2] = 0;
    sector_[dim] = is_hi ? 1 : -1;
  }

  template<typename It>
  void copy(It first) const
  { for (int d=0; d<3; ++d, first++) *first = sector_[d]; }

  bool operator==(const sector_t & other) const
  {
    return
      sector_[0] == other.sector_[0] &&
      sector_[1] == other.sector_[1] &&
      sector_[2] == other.sector_[2];
  }

  int_t nnz() const {
    return
      static_cast<bool>(sector_[0]) +
      static_cast<bool>(sector_[1]) +
      static_cast<bool>(sector_[2]);
  }
  
  void flip(int_t i) { sector_[i] *= -1; }

  void flip() {
    sector_[0] *= -1;
    sector_[1] *= -1;
    sector_[2] *= -1;
  }

  int_t derives_from(const sector_t & other) const
  { 
    return
      std::max(sector_[0]*other[0], 0) +
      std::max(sector_[1]*other[1], 0) +
      std::max(sector_[2]*other[2], 0);
  }

  int_t & operator[](int_t d)
  { return sector_[d]; }
  const int_t & operator[](int_t d) const
  { return sector_[d]; }

  const int_t * begin() const { return &sector_[0]; }
  int_t * begin() { return &sector_[0]; }
  
  const int_t * end() const { return &sector_[3]; }
  int_t * end() { return &sector_[3]; }

  int_t first_non_zero() const
  {
    int_t pos=0;
    for (; pos<3; ++pos)
      if (sector_[pos] != 0)
        break;
    return pos;
  }

  bool is_hi(int_t dim) const
  { return sector_[dim] > 0; }

  friend std::ostream &operator<<( std::ostream &out, const sector_t & sec )
  {
    out << "("  << sec.sector_[0]
        << ", " << sec.sector_[1]
        << ", " << sec.sector_[2] << ")"; 
    return out;
  }

};

////////////////////////////////////////////////////////////////////////////////
/// Flip a sector
////////////////////////////////////////////////////////////////////////////////
inline sector_t flip(sector_t sec)
{
  sec.flip();
  return sec;
}

////////////////////////////////////////////////////////////////////////////////
/// Sector hash function
////////////////////////////////////////////////////////////////////////////////
struct sector_hash_t
{
  std::size_t operator()(const sector_t& k) const
  { return morton_number(3, k[0]+1, k[1]+1, k[2]+1); }
};

////////////////////////////////////////////////////////////////////////////////
/// Sector map
////////////////////////////////////////////////////////////////////////////////
template<typename T>
using sector_map_t = std::unordered_map<sector_t, T, sector_hash_t>;


} // namespace

#endif
