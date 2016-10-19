#pragma once

#include <array>
#include "utils/Globals.hpp"
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

  FindClosestEdge ( const Eigen::VectorXd& searchPoint );

  /**
   * @brief Finds the closest Edge object contained in the given Mesh object.
   *
   * @return True, if a closest Edge object could be found.
   */
  template<typename CONTAINER_T>
  bool operator() ( CONTAINER_T& container );

  const Eigen::VectorXd& getSearchPoint() const;

  bool hasFound() const;

  double getEuclidianDistance();

  mesh::Edge& getClosestEdge();

  const Eigen::VectorXd& getVectorToProjectionPoint() const;

  /// Returns parametric description value (index 0,1) of projected point.
  double getProjectionPointParameter ( int index ) const;

private:

  static logging::Logger _log;

  Eigen::VectorXd _searchPoint;

  double _shortestDistance;

  Eigen::VectorXd _vectorToProjectionPoint;

  std::array<double,2> _parametersProjectionPoint;

  mesh::Edge* _closestEdge;

  void find ( mesh::Edge& edge );
};



// --------------------------------------------------------- HEADER DEFINITIONS

template<typename CONTAINER_T>
bool FindClosestEdge:: operator()
(
  CONTAINER_T& container )
{
  size_t index = 0;
  for ( mesh::Edge& edge : container.edges() ) {
    find ( edge );
    index ++;
  }
  return _closestEdge != nullptr;
}

}} // namespace precice, query


