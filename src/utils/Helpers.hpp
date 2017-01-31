#pragma once

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <set>

namespace precice {
namespace utils {

/// Returns true, if numerical truncation happens in case of type conversion.
template <class Out, class In>
inline bool isTruncated(In in) {
  return (                  in > std::numeric_limits<Out>::max()) ||
         (In(-1) < In(0) && in < std::numeric_limits<Out>::min());
}

/// Exclusive "or" logical operation. Returns true, if either lhs or rhs are true.
inline bool xOR ( bool lhs, bool rhs )
{
   return (lhs && (!rhs)) || ((!lhs) && rhs);
}

/// Evaluates a string to find out if it represents a bool.
/**
 * Returns True if string is yes, true, 1 or on. Otherwise False.
 * This function is case-insensitive.
 */
bool convertStringToBool(std::string const & value);

/**
 * @brief Returns true, if given element is in vector, otherwise false.
 *
 * Requirements:
 * - ELEMENT_T must be comparable by ==
 */
template <typename ELEMENT_T>
bool contained(const ELEMENT_T& element, const std::vector<ELEMENT_T>& vec)
{
  return std::find(vec.begin(), vec.end(), element) != vec.end();
}

template<typename KEY_T, typename ELEMENT_T>
bool contained (
  const KEY_T&                     key,
  const std::map<KEY_T,ELEMENT_T>& map )
{
  return map.find(key) != map.end();
}

template<typename KEY_T>
bool contained (
  const KEY_T&           key,
  const std::set<KEY_T>& set )
{
  return set.find(key) != set.end();
}

/// Returns true if machine is big-endian needed for parallel vtk output
bool isMachineBigEndian();

}} // namespace precice, utils

