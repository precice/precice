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
 * \tparam Compare type of the comparator.
 * \param c the container to check for unique elements.
 * \param comp the comparator to use for the elements.
 * \returns weather all elements in c are unique.
 */
template <typename Container, typename Compare>
bool unique_elements(const Container& c, Compare comp)
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
                [&comp, cbegin](const typename Container::value_type &v) -> bool {
                    return comp(*cbegin, v); 
                }) != cend) {
      return false;
    }
  }
  return true;
}

/** Checks weather the given elements contains no duplicates.
 *
 * \tparam T type of items to check
 * \param c the elements to compare
 * \returns weather all elements are unique.
 */
template <typename Container>
bool unique_elements(const Container& c)
{
  using Comp = std::equal_to<typename Container::value_type>;
  return unique_elements<Container, Comp>(c, Comp{});
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
