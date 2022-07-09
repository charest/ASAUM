#ifndef MERGE_SORT_HPP
#define MERGE_SORT_HPP

namespace prl {

////////////////////////////////////////////////////////////////////////////////
/// merge sort algorithm
////////////////////////////////////////////////////////////////////////////////

namespace detail {

//==============================================================================
// Merge the two sorted arrays together
//==============================================================================
template<typename RandomIt, typename Compare>
void _merge(
  RandomIt d_first,
  RandomIt middle,
  RandomIt d_last,
  Compare comp)
{
  
	using value_type = std::decay_t< decltype(*d_first) >;

  auto pivot = std::next(middle);
  auto end = std::next(d_last);

  // Create temp arrays and copy data to temp arrays L[] and R[]
	std::vector<value_type> L(d_first, pivot), R(pivot, end);

  // Merge the temp arrays back into arr[l..r]

  // Initial index of first subarray
  auto i = L.begin();

  // Initial index of second subarray
  auto j = R.begin();

  // Initial index of merged subarray
  auto k = d_first;

  while (i < L.end() && j < R.end()) {
      if (comp(*j, *i)) {
        *k = *j;
        ++j;
      }
      else {
        *k = *i;
        ++i;
      }
      ++k;
  }

  // Copy the remaining elements of
  // L[], if there are any
  while (i < L.end()) {
      *k = *i;
      ++i;
      ++k;
  }

  // Copy the remaining elements of
  // R[], if there are any
  while (j < R.end()) {
      *k = *j;
      ++j;
      ++k;
  }
}


//==============================================================================
// Main quick sort routine
//==============================================================================
template<typename RandomIt, typename Compare>
void _merge_sort(
  RandomIt d_first,
  RandomIt d_last,
  Compare comp)
{
  // If the subarray is bigger than one element in size...
  if(d_first < d_last) {
    // Partition the array.
    auto n = std::distance(d_first, d_last);
    auto pivot = std::next(d_first, n/2);

    // Recursively divide and partition the two subarrays.
    _merge_sort(d_first, pivot, comp);
    _merge_sort(pivot+1, d_last, comp);
    _merge(d_first, pivot, d_last, comp);
  }
}

} // namespace

//==============================================================================
// Merge sort entry point
//==============================================================================
template<typename RandomIt, typename Compare>
void merge_sort(RandomIt first, RandomIt last, Compare comp)
{
  if (first == last) return;
  detail::_merge_sort(first, std::prev(last), comp);
}

template<typename RandomIt>
void merge_sort(RandomIt first, RandomIt last)
{
  using value_type = std::decay_t< decltype(*first) >;
  merge_sort(first, last, std::less<value_type>());
}

} // namespace

#endif
