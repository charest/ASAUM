#ifndef REORDER_HPP
#define REORDER_HPP

#include <iterator>
#include <utility>

namespace prl {

//! \brief Reorders a vector
template<typename T, typename RandomIt>
void reorder(std::vector<T> & data, RandomIt order)
{
  std::vector<T> new_data(data.size());
  for (size_t i=0; i<data.size(); ++i, ++order)
    new_data[i] = data[ *order ];
  std::swap( new_data, data );
}

//! \brief Reorders an array in place
//! \remark this version maintains the order array
//! \param [in] order_begin The begin iterator for the order array
//! \param [in] order_end   The end iterator for the order array
//! \param [in,out] v The begin iterator for the value array
template<typename order_iterator, typename value_iterator>
void
reorder(const order_iterator order_begin,
  const order_iterator order_end,
  const value_iterator v) {
  using index_t = typename std::iterator_traits<order_iterator>::value_type;
  using diff_t = typename std::iterator_traits<order_iterator>::difference_type;

  diff_t remaining = order_end - 1 - order_begin;
  for(index_t s = index_t(), d; remaining > 0; ++s) {
    if ( s == order_begin[s] ) continue;
    for(d = order_begin[s]; d > s; d = order_begin[d]);
    if(d == s) {
      --remaining;
      auto temp = v[s];
      while(d = order_begin[d], d != s) {
        swap(temp, v[d]);
        --remaining;
      }
      v[s] = temp;
    }
  }
}

//! \brief Reorders an array in place
//! \remark this version destroys the order array for performance gains
//! \param [in,out] order_begin The begin iterator for the order array
//! \param [in,out] order_end   The end iterator for the order array
//! \param [in,out] v The begin iterator for the value array
template<typename order_iterator, typename value_iterator>
void
reorder_destructive(const order_iterator order_begin,
  const order_iterator order_end,
  const value_iterator v) {
  using index_t = typename std::iterator_traits<order_iterator>::value_type;
  using diff_t = typename std::iterator_traits<order_iterator>::difference_type;

  diff_t remaining = order_end - 1 - order_begin;
  for(auto s = index_t(); remaining > 0; ++s) {
    auto d = order_begin[s];
    if(d == diff_t(-1))
      continue;
    --remaining;
    auto temp = v[s];
    for(index_t d2; d != s; d = d2) {
      std::swap(temp, v[d]);
      std::swap(order_begin[d], d2 = index_t(-1));
      --remaining;
    }
    v[s] = temp;
  }
}

} // namespace

#endif
