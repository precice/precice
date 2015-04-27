// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_HELPERS_HPP_
#define PRECICE_UTILS_HELPERS_HPP_

#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <stack>
#include <set>
#include <list>
#include "Parallel.hpp"
#include "Tracer.hpp"

namespace precice {
namespace utils {

/**
 * @brief Returns true, if numerical truncation happens in case of type
 *        conversion.
 */
template <class Out, class In>
inline bool isTruncated(In in) {
  return (                  in > std::numeric_limits<Out>::max()) ||
         (In(-1) < In(0) && in < std::numeric_limits<Out>::min());
}

/**
 * @brief Exclusive "or" logical operation.
 *
 * @brief Returns true, if either lhs or rhs are true.
 */
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
const bool contained(const ELEMENT_T& element, const std::vector<ELEMENT_T>& vec)
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

/**
 * @brief Enables overload of operator() and operator, to fill STL containers.
 */
template<typename CONTAINER_T>
struct AppendIterator
{
public:

  /**
   * @brief Constructor, used by appedTo and operator+.
   */
  AppendIterator ( CONTAINER_T& container )
  :
    _container ( container )
  {};

  /**
   * @brief Enables to chain assignments such as (1)(2)(3).
   */
  template<typename CONTAINED_T>
  AppendIterator& operator() ( const CONTAINED_T& right )
  {
    _container.insert ( _container.end(), right );
    return *this;
  }

  /**
   * @brief Enables to chain assignments such as 1, 2, 3.
   */
  template<typename CONTAINED_T>
  AppendIterator& operator, ( const CONTAINED_T& right )
  {
    _container.insert ( _container.end(), right );
    return *this;
  }

private:

  // @brief Container class assigned to.
  CONTAINER_T& _container;
};

/**
 * @brief Enables to append values to a STL conform container.
 *
 * ! Examples
 * std::list<double> _list;
 * appendTo(_list) (1.0)(2.0)(3.0);
 * appendTo(_list) 4.0, 5.0, 6.0;
 *
 * ! Requirements
 * - CONTAINER_T must offer a function insert(iterator_t i,
 *                                            const CONTAINED_T & item)
 */
template<typename CONTAINER_T>
inline AppendIterator<CONTAINER_T> appendTo
(
   CONTAINER_T& left )
{
   return AppendIterator<CONTAINER_T>(left);
}

template <class T, class Container = std::deque<T> >
class Stack : public std::stack<T, Container> {
public:
  typedef typename Container::iterator iterator;
  typedef typename Container::const_iterator const_iterator;

  iterator
  begin() {
    return this->c.begin();
  }

  iterator
  end() {
    return this->c.end();
  }

  const_iterator
  begin() const {
    return this->c.begin();
  }

  const_iterator
  end() const {
    return this->c.end();
  }
};

}} // namespace precice, utils

// ------------------------------------------------------------- FREE FUNCTIONS

namespace precice {

/**
 * @brief Enables to append values to a std::vector.
 *
 * ! Examples
 * <code>
 * std::vector<int> _vector;
 * _vector += 1, 2, 3;
 * _vector += (1)(2)(3);
 * </code>
 *
 * ! Requirements
 * - CONTAINER_T must offer a function insert(iterator_t i,
 *                                            const CONTAINED_T & item)
 */
template< typename CONTAINED_T >
inline utils::AppendIterator< std::vector<CONTAINED_T> > operator+=
(
   std::vector<CONTAINED_T> & left,
   const CONTAINED_T        & right )
{
   left.insert ( left.end(), right );
   return utils::AppendIterator< std::vector<CONTAINED_T> > ( left );
}

} // namespace precice

#endif /* PRECICE_UTILS_HELPERS_HPP_ */
