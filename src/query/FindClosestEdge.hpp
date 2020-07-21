#pragma once

#include <Eigen/Core>
#include <array>
#include <limits>
#include <stddef.h>
#include <string>
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
class Edge;
class Mesh;
} // namespace mesh
} // namespace precice

namespace precice {
namespace query {

/// Finds the closest Edge object contained in Mesh objects.
class FindClosestEdge {
public:
  explicit FindClosestEdge(const Eigen::VectorXd &searchPoint);

  /**
   * @brief Finds the closest Edge object contained in the given Mesh object.
   *
   * @return True, if a closest Edge object could be found.
   */
  template <typename CONTAINER_T>
  bool operator()(CONTAINER_T &container);

  const Eigen::VectorXd &getSearchPoint() const;

  bool hasFound() const;

  double getEuclidianDistance();

  mesh::Edge &getClosestEdge();

  const Eigen::VectorXd &getVectorToProjectionPoint() const;

  /// Returns parametric description value (index 0,1) of projected point.
  double getProjectionPointParameter(int index) const;

private:
  logging::Logger _log{"query::FindClosestEdge"};

  Eigen::VectorXd _searchPoint;

  double _shortestDistance = std::numeric_limits<double>::max();

  Eigen::VectorXd _vectorToProjectionPoint;

  std::array<double, 2> _parametersProjectionPoint;

  mesh::Edge *_closestEdge = nullptr;

  void find(mesh::Edge &edge);
};

// --------------------------------------------------------- HEADER DEFINITIONS

template <typename CONTAINER_T>
bool FindClosestEdge::operator()(
    CONTAINER_T &container)
{
  size_t index = 0;
  for (mesh::Edge &edge : container.edges()) {
    find(edge);
    index++;
  }
  return _closestEdge != nullptr;
}

} // namespace query
} // namespace precice
