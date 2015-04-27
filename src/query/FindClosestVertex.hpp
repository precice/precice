// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QUERY_FINDCLOSESTVERTEX_HPP_
#define PRECICE_QUERY_FINDCLOSESTVERTEX_HPP_

#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "mesh/Vertex.hpp"
#include <limits>

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/**
 * @brief Finds the closest Vertex object hold by Mesh objects.
 */
class FindClosestVertex
{
public:

  /**
   * @brief Constructor.
   *
   * @param searchPoint [IN] Coordinates of origin of search for closest point.
   */
  template<typename VECTOR_T>
  FindClosestVertex ( const VECTOR_T& searchPoint );

  /**
   * @brief Searches among all Vertex objects hold by the given Mesh object.
   *
   * When called for different meshes, the closest Vertex object in among all
   * the meshes is found.
   */
  template<typename CONTAINER_T>
  bool operator() ( CONTAINER_T& container );

  /**
   * @brief Returns the coordinates of the search point.
   */
  const utils::DynVector& getSearchPoint() const;

  bool hasFound() const;

  /**
   * @brief Returns the euclidian distance to the closest vertex found.
   *
   * Precondition: find() has been called and returned true.
   */
  double getEuclidianDistance();

  /**
   * @brief Returns the closest vertex found.
   *
   * Precondition: find() has been called and returned true.
   */
  mesh::Vertex& getClosestVertex();

private:

  // @brief Origin of search for closest vertex.
  utils::DynVector _searchPoint;

  // @brief Distance to closest Vertex object found.
  double _shortestDistance;

  // @brief Pointer to closest Vertex object found.
  mesh::Vertex* _closestVertex;
};

// --------------------------------------------------------- HEADER DEFINITIONS

template<typename VECTOR_T>
FindClosestVertex:: FindClosestVertex
(
  const VECTOR_T& searchPoint )
:
  _searchPoint (searchPoint),
  _shortestDistance (std::numeric_limits<double>::max()),
  _closestVertex (NULL)
{}

template<typename CONTAINER_T>
bool FindClosestVertex:: operator()
(
  CONTAINER_T& container )
{
  utils::DynVector vectorDistance(_searchPoint.size(), 0.0);
  foreach ( mesh::Vertex& vertex, container.vertices() ) {
    assertion2 ( vertex.getDimensions() == _searchPoint.size(),
                 vertex.getDimensions(), _searchPoint.size() );
    vectorDistance = vertex.getCoords();
    vectorDistance -= _searchPoint;
    double distance = tarch::la::norm2(vectorDistance);
    if ( distance < _shortestDistance) {
      _shortestDistance = distance;
      _closestVertex = &vertex;
    }
  }
  return _closestVertex != NULL;
}

}} // namespace precice, query

#endif /* PRECICE_QUERY_FINDCLOSESTVERTEX_HPP_ */
