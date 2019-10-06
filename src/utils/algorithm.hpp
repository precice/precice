#pragma once

#include<algorithm>
#include<array>
#include<functional>
#include<type_traits>

namespace precice {
namespace utils {

/// Function that generates an array from given elements.
template<typename... Elements>
auto make_array(Elements&&... elements) -> std::array<typename std::common_type<Elements...>::type, sizeof...(Elements)> {
    return { std::forward<Elements>(elements)... };
}

/** Checks weather the given elements contains no duplicates.
 *
 * \tparam Container type of the passed container.
 * \tparam BinaryPredicate the predicate used to compare two elements for equality.
 * \param c the container to check for unique elements.
 * \returns weather all elements in c are unique.
 */
template <typename Container, typename BinaryPredicate = std::equal_to<typename Container::value_type>>
bool unique_elements(const Container& c, BinaryPredicate p = {})
{
  // This algorithm runs on small containers, so the O(n^2) does not hurt.
  auto cbegin = c.begin();
  auto cend   = c.end();
  // An empty set has unique elements
  if (cbegin == cend) {
    return true;
  }
  auto cstart = cbegin + 1;
  for (; cstart < cend; ++cbegin, ++cstart) {
    if (std::find_if(cstart, cend,
                [&p, cbegin](const typename Container::value_type &v) -> bool {
                    return p(*cbegin, v); 
                }) != cend) {
      return false;
    }
  }
  return true;
}


/** std::mismatch
 * @todo{Remove when migrating to c++14}
 */
template<class InputIt1, class InputIt2>
std::pair<InputIt1, InputIt2>
    mismatch(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2)
{
    while (first1 != last1 && first2 != last2 && *first1 == *first2) {
        ++first1, ++first2;
    }
    return std::make_pair(first1, first2);
}


} // namespace utils
} // namespace precice
