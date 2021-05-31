#pragma once

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace precice {
namespace utils {

/// Returns true, if numerical truncation happens in case of type conversion.
template <class Out, class In>
inline bool isTruncated(In in)
{
  return (in > std::numeric_limits<Out>::max()) ||
         (In(-1) < In(0) && in < std::numeric_limits<Out>::min());
}

/// Returns true if the argument represents a vaild port
inline bool isValidPort(int port)
{
  return (port >= 0) && !utils::isTruncated<unsigned short>(port);
}

/// Exclusive "or" logical operation. Returns true, if either lhs or rhs are true.
inline bool xOR(bool lhs, bool rhs)
{
  return (lhs && (!rhs)) || ((!lhs) && rhs);
}

/**
 * @brief Returns true, if given element is in vector, otherwise false.
 *
 * Requirements:
 * - ELEMENT_T must be comparable by ==
 */
template <typename ELEMENT_T>
bool contained(const ELEMENT_T &element, const std::vector<ELEMENT_T> &vec)
{
  return std::find(vec.begin(), vec.end(), element) != vec.end();
}

template <typename KEY_T, typename ELEMENT_T>
bool contained(
    const KEY_T &                     key,
    const std::map<KEY_T, ELEMENT_T> &map)
{
  return map.find(key) != map.end();
}

template <typename KEY_T>
bool contained(
    const KEY_T &          key,
    const std::set<KEY_T> &set)
{
  return set.find(key) != set.end();
}

/// Returns true if machine is big-endian needed for parallel vtk output
bool isMachineBigEndian();

} // namespace utils
} // namespace precice
