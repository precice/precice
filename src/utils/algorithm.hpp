#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>

namespace precice {
namespace utils {

/// Function that generates an array from given elements.
template <typename... Elements>
auto make_array(Elements &&... elements) -> std::array<typename std::common_type<Elements...>::type, sizeof...(Elements)>
{
  return {std::forward<Elements>(elements)...};
}

/** Checks weather the given elements contains no duplicates.
 *
 * \tparam Container type of the passed container.
 * \tparam BinaryPredicate the predicate used to compare two elements for equality.
 * \param c the container to check for unique elements.
 * \returns weather all elements in c are unique.
 */
template <typename Container, typename BinaryPredicate = std::equal_to<typename Container::value_type>>
bool unique_elements(const Container &c, BinaryPredicate p = {})
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

/** intersperse a the range [first, last[ with a given element.
 * 
 * This results in a range [first, elem, first+1, elem, ... , elem, last[
 * 
 * \tparam InputIter the type of the input iterators
 * \tparam ElemT the type of the element to intersperse
 */
template <class InputIter, class ElemT>
void intersperse(InputIter first, InputIter last, const ElemT &elem, std::ostream &out)
{
  if (first == last)
    return;

  out << *first++;
  for (; first != last; ++first) {
    out << elem << *first;
  }
}

/** std::mismatch
 * @todo{Remove when migrating to c++14}
 */
template <class InputIt1, class InputIt2>
std::pair<InputIt1, InputIt2>
mismatch(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2)
{
  while (first1 != last1 && first2 != last2 && *first1 == *first2) {
    ++first1, ++first2;
  }
  return std::make_pair(first1, first2);
}

/// The RangePreview object used as a lazy proxy struct for proviewing the content of a Range
template <typename InputIter>
struct RangePreview {
  using Size = typename std::iterator_traits<InputIter>::difference_type;
  Size      n;
  InputIter begin;
  InputIter end;

  RangePreview(Size n, InputIter begin, InputIter end)
      : n(n), begin(begin), end(end) {}

  void print(std::ostream &out) const
  {
    if (begin == end) {
      out << "<Empty Range>";
      return;
    }

    out << '[';

    if (n == 0) {
      out << " ... ";
    } else {
      auto              dist = std::distance(begin, end);
      const std::string sep{", "};
      if (dist <= n * 2) {
        intersperse(begin, end, sep, out);
      } else {
        auto last1 = begin;
        std::advance(last1, n);
        intersperse(begin, last1, sep, out);
        out << sep << "... " << sep;
        auto first2 = begin;
        std::advance(first2, dist - n);
        intersperse(first2, end, sep, out);
      }
    }

    auto mm = std::minmax_element(begin, end);
    out << "] min:" << *mm.first << " max:" << *mm.second;
  }
};

/// Allows streaming of RangePreview objects
template <typename Iter>
std::ostream &operator<<(std::ostream &out, const RangePreview<Iter> &rp)
{
  rp.print(out);
  return out;
}

/** returns a display object which previews a range
 *
 * The preview contains the first and last n elements and the minmax-elements.
 */
template <typename Range, typename Iter = typename Range::const_iterator, typename Size = typename std::iterator_traits<Iter>::difference_type>
const RangePreview<Iter> previewRange(Size n, const Range &range)
{
  return {n, std::begin(range), std::end(range)};
}

/// Reorders an array given an array of unique indices.
template <typename T, typename Index, size_t n>
auto reorder_array(const std::array<Index, n> &order, const std::array<T, n> &elements) -> std::array<T, n>
{
  std::array<T, n> reordered;
  for (std::size_t i = 0; i < n; ++i) {
    reordered[i] = elements[order[i]];
  }
  return reordered;
}

} // namespace utils
} // namespace precice
