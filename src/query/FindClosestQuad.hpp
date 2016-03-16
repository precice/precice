// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QUERY_FINDCLOSESTQUAD_HPP_
#define PRECICE_QUERY_FINDCLOSESTQUAD_HPP_

#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include "boost/array.hpp"
#include "boost/assign.hpp"
#include <limits>

namespace precice {
  namespace mesh {
    class Mesh;
    class Quad;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/**
 * @brief Finds the closest Quad object contained in a Mesh object.
 */
class FindClosestQuad
{
public:

  /**
   * @brief Constructor.
   *
   * @param searchPoint [IN] Origin from which the closest Quad object
   *        should be found.
   */
  template<typename VECTOR_T>
  FindClosestQuad ( const VECTOR_T& searchPoint );

  /**
   * @brief Searches for the closest quad on the given mesh object.
   *
   * @return True, if a closest quad could be found.
   */
  template<typename CONTAINER_T>
  bool operator() ( CONTAINER_T & container );

  /**
   * @brief Returns the coordinates of the search point.
   */
  const utils::DynVector& getSearchPoint() const;

  /**
   * @brief Returns true, if a closest quad has been found.
   */
  bool hasFound() const;

  /**
   * @brief Returns the distance to the found quad.
   *
   * Precondition: Find has been called and returned true.
   */
  double getEuclidianDistance();

  /**
   * @brief Returns the found quad.
   *
   * Precondition: Find has been called and returned true.
   */
  mesh::Quad& getClosestQuad();

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

  static logging::Logger _log;

  // @brief Search point coordinates.
  utils::DynVector _searchPoint;

  // @brief Shortest distance to the found Quad object.
  double _shortestDistance;

  // @brief Vector from search point to projection point.
  utils::DynVector _vectorToProjectionPoint;

  // @brief Quad coordinates of the projection point.
  boost::array<double,4> _parametersProjectionPoint; // Does this make sense?

  // @brief Pointer to found Quad object.
  mesh::Quad* _closestQuad;

  void find ( mesh::Quad& quad );
};

// --------------------------------------------------------- HEADER DEFINITIONS

template<typename VECTOR_T>
FindClosestQuad:: FindClosestQuad
(
  const VECTOR_T& searchPoint )
:
  _searchPoint ( searchPoint ),
  _shortestDistance ( std::numeric_limits<double>::max() ),
  _vectorToProjectionPoint ( _searchPoint.size(), std::numeric_limits<double>::max() ),
  _parametersProjectionPoint ( boost::assign::list_of
      (_shortestDistance)(_shortestDistance)(_shortestDistance)(_shortestDistance)
      .to_array(_parametersProjectionPoint)),
  _closestQuad ( NULL )
{}

template<typename CONTAINER_T>
bool FindClosestQuad:: operator() ( CONTAINER_T& container )
{
  foreach ( mesh::Quad& quad, container.quads() ) {
    find ( quad );
  }
  return _closestQuad != NULL;
}

}} // namespace precice, query

#endif /* PRECICE_QUERY_FINDCLOSESTQUAD_HPP_ */
