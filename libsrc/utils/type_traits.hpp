#ifndef TYPE_TRAITS_HPP
#define TYPE_TRAITS_HPP

#include <vector>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
//! \brief Check to see if a type is a vector
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct is_vector {
  static constexpr bool value = false;
};

template <typename T>
struct is_vector<std::vector<T>> {
  static constexpr bool value = true;
};

////////////////////////////////////////////////////////////////////////////////
//! \brief Check if a particular type T is a container.
//! \remark If T is not, this version is instantiated.
//! \remark This version adheres to the strict requirements of an STL container.
////////////////////////////////////////////////////////////////////////////////
namespace detail {
template <typename... Ts>
struct is_container {
};
} // namespace

template <typename T, typename _ = void>
struct is_container : std::false_type {
};

//! \brief Check if a particular type T is a container.
//! \remark If T is, this version is instantiated.
//! \remark This version adheres to the strict requirements of an STL container.
template <typename T>
struct is_container<T,
  std::conditional_t<false,
    detail::is_container<typename T::value_type, typename T::size_type,
      typename T::iterator,
      typename T::const_iterator, decltype(std::declval<T>().size()),
      decltype(std::declval<T>().begin()), decltype(std::declval<T>().end()),
      decltype(std::declval<T>().cbegin()), decltype(std::declval<T>().cend()),
      decltype(std::declval<T>().data())>,
    void>> : public std::true_type {
};

//! \brief Equal to true if T is a container.
//! \remark This version adheres to the strict requirements of an STL container.
template <typename T>
constexpr bool is_container_v = is_container<T>::value;

} // namespace

#endif
