#ifndef INSERTION_SORT_HPP
#define INSERTION_SORT_HPP

namespace prl {

////////////////////////////////////////////////////////////////////////////////
///  sort algorithm
////////////////////////////////////////////////////////////////////////////////

namespace detail {

//==============================================================================
// Main insertion sort routine
//==============================================================================
template<typename RandomIt, typename Compare>
void _insert_sort(
  RandomIt d_first,
  RandomIt d_last,
  Compare comp)
{
  // If the array is bigger than one element in size...
  if (d_first < d_last) {

    auto end = std::next(d_last);
    RandomIt j;
    for (auto i=std::next(d_first); i < end; ++i)
    {
      j = std::prev(i);
      auto key = *i;

      // Move elements of arr[0..i-1], that are
      // greater than key, to one position ahead
      // of their current position
      while (j >= d_first && comp(key, *j))
      {
        std::swap(*j, *std::next(j));
        --j;
      }

      // swap
      *std::next(j) = key;
    }

  } // no-op
}

} // namespace

//==============================================================================
// Insert sort entry point
//==============================================================================
template<typename RandomIt, typename Compare>
void insert_sort(RandomIt first, RandomIt last, Compare comp)
{
  if (first == last) return;
  detail::_insert_sort(first, std::prev(last), comp);
}

template<typename RandomIt>
void insert_sort(RandomIt first, RandomIt last)
{
  using value_type = std::decay_t< decltype(*first) >;
  insert_sort(first, last, std::less<value_type>());
}

} // namespace

#endif
