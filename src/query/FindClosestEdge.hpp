// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QUERY_FINDCLOSESTEDGE_HPP_
#define PRECICE_QUERY_FINDCLOSESTEDGE_HPP_

#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "boost/array.hpp"
#include "utils/Globals.hpp"
#include "boost/assign.hpp"
#include <limits>

namespace precice {
  namespace mesh {
    class Edge;
    class Mesh;
  }
}

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/**
 * @brief Finds the closest Edge object contained in Mesh objects.
 */
class FindClosestEdge
{
public:

  template<typename VECTOR_T>
  FindClosestEdge ( const VECTOR_T& searchPoint );

  /**
   * @brief Finds the closest Edge object contained in the given Mesh object.
   *
   * @return True, if a closest Edge object could be found.
   */
  template<typename CONTAINER_T>
  bool operator() ( CONTAINER_T& container );

  const utils::DynVector& getSearchPoint() const;

  bool hasFound() const;

  double getEuclidianDistance();

  mesh::Edge& getClosestEdge();

  const utils::DynVector& getVectorToProjectionPoint() const;

  /**
   * @brief Returns parametric description value (index 0,1) of projected point.
   */
  double getProjectionPointParameter ( int index ) const;

private:

  static tarch::logging::Log _log;

  utils::DynVector _searchPoint;

  double _shortestDistance;

  utils::DynVector _vectorToProjectionPoint;

  boost::array<double,2> _parametersProjectionPoint;

  mesh::Edge* _closestEdge;

  void find ( mesh::Edge& edge );
};



// --------------------------------------------------------- HEADER DEFINITIONS

template<typename VECTOR_T>
FindClosestEdge:: FindClosestEdge
(
  const VECTOR_T& searchPoint )
:
  _searchPoint ( searchPoint ),
  _shortestDistance ( std::numeric_limits<double>::max() ),
  _vectorToProjectionPoint ( searchPoint.size(), std::numeric_limits<double>::max() ),
  _parametersProjectionPoint (
      boost::assign::list_of(_shortestDistance)(_shortestDistance).to_array(
      _parametersProjectionPoint)),
  _closestEdge ( NULL )
{
  assertion1 ( (_searchPoint.size() == 2) || (_searchPoint.size() == 3),
               _searchPoint.size() );
}

template<typename CONTAINER_T>
bool FindClosestEdge:: operator()
(
  CONTAINER_T& container )
{
  size_t index = 0;
  foreach ( mesh::Edge& edge, container.edges() ) {
    find ( edge );
    index ++;
  }
  return _closestEdge != NULL;
}

}} // namespace precice, query

#endif /* PRECICE_QUERY_FINDCLOSESTEDGE_HPP_ */
