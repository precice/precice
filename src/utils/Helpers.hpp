#pragma once

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include "Parallel.hpp"

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

///**
// * @brief Returns the corresponding zero value/object to value_t.
// */
//template< typename value_t >
//inline value_t getZero() { return value_t(0.0); }
//
///**
// * @brief Template specialization to handle primitive type double.
// */
//template<>
//inline double getZero<double>() { return 0.0; };
//
///**
// * @brief Template specialization to handle primitive type int.
// */
//template<>
//inline int getZero<int>() { return 0; };

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

std::string getTypeName(const double& var);

std::string getTypeName(const std::string& var);

std::string getTypeName(const bool& var);

std::string getTypeName(const int& var);

/// Returns true if machine is big-endian needed for parallel vtk output
bool isMachineBigEndian();

}} // namespace precice, utils

