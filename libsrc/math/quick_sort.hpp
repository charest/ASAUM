#ifndef QUICK_SORT_HPP
#define QUICK_SORT_HPP

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// quick sort algorithm
////////////////////////////////////////////////////////////////////////////////

namespace detail {

//==============================================================================
// Partition algorithm
//==============================================================================
template<typename RandomIt, typename Compare>
RandomIt _partition(
  RandomIt d_first,
  RandomIt d_last,
  Compare comp)
{
  // Choose our pivot as the last element in the array.
  auto pivot = *d_last;
  
	// index of smaller element
  auto i = d_first-1;
  
  // Go through each element in this subarray.
  for(auto j = d_first; j < d_last; ++j) {

    // If the current element is smaller than (or equal to) the pivot...
    if(!comp(pivot, *j)) {
      // Increment index of the smaller element
      ++i;
      // swap
      std::swap(*i, *j);
    }

  }
  
  // final swap
  i++;
  std::swap(*i, *d_last);
  
  // Return the index of the pivot.
  return i;
}


//==============================================================================
// Main quick sort routine
//==============================================================================
template<typename RandomIt, typename Compare>
void _quick_sort(
  RandomIt d_first,
  RandomIt d_last,
  Compare comp)
{
  // If the subarray is bigger than one element in size...
  if(d_first < d_last) {
    // Partition the array.
    auto pivot = _partition(d_first, d_last, comp);
      
    // Recursively divide and partition the two subarrays.
    _quick_sort(d_first, pivot-1, comp);
    _quick_sort(pivot+1, d_last, comp);
  }
}

} // namespace

//==============================================================================
// Quick sort entry point
//==============================================================================
template<typename RandomIt, typename Compare>
void quick_sort(RandomIt first, RandomIt last, Compare comp)
{
  if (first == last) return;
  detail::_quick_sort(first, std::prev(last), comp);
}

template<typename RandomIt>
void quick_sort(RandomIt first, RandomIt last)
{
  using value_type = std::decay_t< decltype(*first) >;
  quick_sort(first, last, std::less<value_type>());
}

} // namespace

#endif
