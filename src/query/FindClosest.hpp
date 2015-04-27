// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QUERY_FINDCLOSEST_HPP_
#define PRECICE_QUERY_FINDCLOSEST_HPP_

#include "FindClosestVertex.hpp"
#include "FindClosestEdge.hpp"
#include "FindClosestTriangle.hpp"
#include "FindClosestQuad.hpp"
#include "mesh/PropertyContainer.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"
#include <map>

namespace precice {
   namespace mesh {
      class Mesh;
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/**
 * @brief Weighting and reference to target element for a value to interpolate
 */
struct InterpolationElement
{
  mesh::Vertex* element;
  double weight;

  InterpolationElement ()
  : element(NULL), weight(0.0) {}
};

/**
 * @brief Closest element to all objects with given geometry ID
 */
struct ClosestElement
{
  std::vector<int> meshIDs;
  double distance;
  utils::DynVector vectorToElement;
  std::vector<InterpolationElement> interpolationElements;

//  ClosestElement ()
//  : meshIDs(), distance(0.0), vectorToElement(), interpolationElements() {}

  ClosestElement (int dim)
  : meshIDs(), distance(0.0), vectorToElement(dim,0.0), interpolationElements() {}
};


/**
 * @brief Determines closest Triangle, Edge, or Vertex object to a given point.
 *
 * Computes a distance vector to every Triangle and Vertex object found and
 * stores the object with shortest distance.
 *
 * The distance to a Vertex object is measured directly from the search point
 * to the Vertex's coordinates. The interpolation weight to a vertex is one.
 *
 * The distance to a Triangle is measured from the search point, to a point
 * orthogonally projected onto the triangles plane (or line in 2D). The
 * distance is valid only then, when the projected point lies within the
 * triangle, i.e. when the barycentric coordinates of the projected point are
 * all smaller or equal to one. The interpolation weigths for the points of the
 * triangle are equal to the barycentric coordinates.
 */
class FindClosest
{
public:

  /**
   * @brief Constructor, searchpoint can be specified only there
   *
   * @param searchpoint [IN] Point from where distances to objects are measured
   */
  template<typename VECTOR_T>
  FindClosest ( const VECTOR_T& searchpoint );

  /// Finds closest distance to all mesh elements in the given container.
  template<typename CONTAINER_T>
  bool operator() ( CONTAINER_T& container );

  /// Returns true, if a closest element was found.
  bool hasFound() const;

  /// Returns ClosestElement found, error when no visitable has been found
  const ClosestElement& getClosest();

  /// Returns the euclidian distance to the closest element.
  double getEuclidianDistance();

  /// Returns search point
  const utils::DynVector& getSearchPoint() const;

  /// Resets the found visitables, not done automatically
  void reset();

private:

  static tarch::logging::Log _log;

  /// Finds closest distance to Vertex objects.
  FindClosestVertex _findClosestVertex;

  /// Finds closest distance to Edge objects.
  FindClosestEdge _findClosestEdge;

  /// Find closest distance to Triangle objects.
  FindClosestTriangle _findClosestTriangle;

  /// Find closest distance to Quad objects.
  FindClosestQuad _findClosestQuad;

  /// Closest mesh element.
  ClosestElement _closest;

  /// Search point, from where distances to objects are measured
  const utils::DynVector _searchpoint;

  /**
   * @brief Determines the closest element from all FindXY member objects.
   *
   * @return True, if a closest object has been found.
   */
  bool determineClosest();
};

// --------------------------------------------------------- HEADER DEFINITIONS

template<typename VECTOR_T>
FindClosest:: FindClosest
(
  const VECTOR_T& searchpoint )
:
  _findClosestVertex(searchpoint),
  _findClosestEdge(searchpoint),
  _findClosestTriangle(searchpoint),
  _findClosestQuad(searchpoint),
  _closest(searchpoint.size()),
  _searchpoint(searchpoint)
{}

template<typename CONTAINER_T>
bool FindClosest:: operator()
(
  CONTAINER_T& container )
{
  // It is not valid here, to stop the search for vertices (e.g.) if a closest
  // edge has been found already, since the edge might be oriented wrongly.
  _findClosestTriangle(container);
  _findClosestEdge(container);
  _findClosestVertex(container);
  _findClosestQuad(container);
  return determineClosest();
}

}} // namespace precice, query

#endif /* PRECICE_QUERY_FINDCLOSEST_HPP_ */
