// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QUERY_FINDCLOSESTTRIANGLE_HPP_
#define PRECICE_QUERY_FINDCLOSESTTRIANGLE_HPP_

#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include "boost/array.hpp"
#include "boost/assign.hpp"
#include <limits>

namespace precice {
  namespace mesh {
    class Mesh;
    class Triangle;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/**
 * @brief Finds the closest Triangle object contained in a Mesh object.
 */
class FindClosestTriangle
{
public:

  /**
   * @brief Constructor.
   *
   * @param searchPoint [IN] Origin from which the closest Triangle object
   *        should be found.
   */
  template<typename VECTOR_T>
  FindClosestTriangle ( const VECTOR_T& searchPoint );

  /**
   * @brief Searches for the closest triangle on the given mesh object.
   *
   * @return True, if a closest triangle could be found.
   */
  template<typename CONTAINER_T>
  bool operator() ( CONTAINER_T & container );

  /**
   * @brief Returns the coordinates of the search point.
   */
  const utils::DynVector& getSearchPoint() const;

  /**
   * @brief Returns true, if a closest triangle has been found.
   */
  bool hasFound() const;

  /**
   * @brief Returns the distance to the found triangle.
   *
   * Precondition: Find has been called and returned true.
   */
  double getEuclidianDistance();

  /**
   * @brief Returns the found triangle.
   *
   * Precondition: Find has been called and returned true.
   */
  mesh::Triangle& getClosestTriangle();

  /**
   * @brief Returns the vector from the search point to the projection point.
   *
   * Precondition: Find has been called and returned true.
   */
  const utils::DynVector& getVectorToProjectionPoint() const;

  /**
   * @brief Returns parametric description value (index 0, 1, 2) of proj. point.
   */
  double getProjectionPointParameter ( int index ) const;

private:

  static tarch::logging::Log _log;

  // @brief Search point coordinates.
  utils::DynVector _searchPoint;

  // @brief Shortest distance to the found Triangle object.
  double _shortestDistance;

  // @brief Vector from search point to projection point.
  utils::DynVector _vectorToProjectionPoint;

  // @brief Barycentric coordinates of the projection point.
  boost::array<double,3> _parametersProjectionPoint;

  // @brief Pointer to found Triangle object.
  mesh::Triangle* _closestTriangle;

  void find ( mesh::Triangle& triangle );
};

// --------------------------------------------------------- HEADER DEFINITIONS

template<typename VECTOR_T>
FindClosestTriangle:: FindClosestTriangle
(
  const VECTOR_T& searchPoint )
:
  _searchPoint ( searchPoint ),
  _shortestDistance ( std::numeric_limits<double>::max() ),
  _vectorToProjectionPoint ( _searchPoint.size(), std::numeric_limits<double>::max() ),
  _parametersProjectionPoint ( boost::assign::list_of(_shortestDistance)
      (_shortestDistance)(_shortestDistance).to_array(_parametersProjectionPoint)),
  _closestTriangle ( NULL )
{}

template<typename CONTAINER_T>
bool FindClosestTriangle:: operator() ( CONTAINER_T& container )
{
  foreach ( mesh::Triangle& triangle, container.triangles() ) {
    find ( triangle );
  }
  return _closestTriangle != NULL;
}

}} // namespace precice, query

#endif /* PRECICE_QUERY_FINDCLOSESTTRIANGLE_HPP_ */
