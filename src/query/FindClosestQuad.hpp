#pragma once

#include <Eigen/Core>
#include <array>
#include <limits>
#include <stddef.h>
#include <string>
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
class Mesh;
class Quad;
} // namespace mesh
} // namespace precice

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/**
 * @brief Finds the closest Quad object contained in a Mesh object.
 */
class FindClosestQuad {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] searchPoint Origin from which the closest Quad object should be found.
   */
  explicit FindClosestQuad(const Eigen::VectorXd &searchPoint);

  /**
   * @brief Searches for the closest quad on the given mesh object.
   *
   * @return True, if a closest quad could be found.
   */
  template <typename CONTAINER_T>
  bool operator()(CONTAINER_T &container);

  /// Returns the coordinates of the search point.
  const Eigen::VectorXd &getSearchPoint() const;

  /// Returns true, if a closest quad has been found.
  bool hasFound() const;

  /**
   * @brief Returns the distance to the found quad.
   *
   * @pre Find has been called and returned true.
   */
  double getEuclidianDistance();

  /**
   * @brief Returns the found quad.
   *
   * @pre Find has been called and returned true.
   */
  mesh::Quad &getClosestQuad();

  /**
   * @brief Returns the vector from the search point to the projection point.
   *
   * @pre Find has been called and returned true.
   */
  const Eigen::VectorXd &getVectorToProjectionPoint() const;

  /// Returns parametric description value (index 0, 1, 2) of proj. point.
  double getProjectionPointParameter(int index) const;

private:
  logging::Logger _log{"query::FindClosestQuad"};

  /// Search point coordinates.
  Eigen::VectorXd _searchPoint;

  /// Shortest distance to the found Quad object.
  double _shortestDistance = std::numeric_limits<double>::max();

  /// Vector from search point to projection point.
  Eigen::VectorXd _vectorToProjectionPoint;

  /// Quad coordinates of the projection point.
  std::array<double, 4> _parametersProjectionPoint; // Does this make sense?

  /// Pointer to found Quad object.
  mesh::Quad *_closestQuad = nullptr;

  void find(mesh::Quad &quad);
};

// --------------------------------------------------------- HEADER DEFINITIONS

template <typename CONTAINER_T>
bool FindClosestQuad::operator()(CONTAINER_T &container)
{
  for (mesh::Quad &quad : container.quads()) {
    find(quad);
  }
  return _closestQuad != NULL;
}

} // namespace query
} // namespace precice
