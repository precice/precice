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
class Triangle;
} // namespace mesh
} // namespace precice

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/// Finds the closest Triangle object contained in a Mesh object.
class FindClosestTriangle {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] searchPoint Origin from which the closest Triangle object should be found.
   */
  explicit FindClosestTriangle(const Eigen::VectorXd &searchPoint);

  /**
   * @brief Searches for the closest triangle on the given mesh object.
   *
   * @return True, if a closest triangle could be found.
   */
  template <typename CONTAINER_T>
  bool operator()(CONTAINER_T &container);

  /// Returns the coordinates of the search point.
  const Eigen::VectorXd &getSearchPoint() const;

  /// Returns true, if a closest triangle has been found.
  bool hasFound() const;

  /**
   * @brief Returns the distance to the found triangle.
   *
   * @pre Find has been called and returned true.
   */
  double getEuclidianDistance() const;

  /**
   * @brief Returns the found triangle.
   *
   * @pre Find has been called and returned true.
   */
  mesh::Triangle &getClosestTriangle();

  /**
   * @brief Returns the vector from the search point to the projection point.
   *
   * @pre Find has been called and returned true.
   */
  const Eigen::VectorXd &getVectorToProjectionPoint() const;

  /// Returns parametric description value (index 0, 1, 2) of proj. point.
  double getProjectionPointParameter(int index) const;

private:
  logging::Logger _log{"query::FindClosestTriangle"};

  /// Search point coordinates.
  Eigen::VectorXd _searchPoint;

  /// Shortest distance to the found Triangle object.
  double _shortestDistance = std::numeric_limits<double>::max();

  /// Vector from search point to projection point.
  Eigen::VectorXd _vectorToProjectionPoint;

  /// Barycentric coordinates of the projection point.
  std::array<double, 3> _parametersProjectionPoint;

  /// Pointer to found Triangle object.
  mesh::Triangle *_closestTriangle = nullptr;

  void find(mesh::Triangle &triangle);
};

// --------------------------------------------------------- HEADER DEFINITIONS

template <typename CONTAINER_T>
bool FindClosestTriangle::operator()(CONTAINER_T &container)
{
  for (mesh::Triangle &triangle : container.triangles()) {
    find(triangle);
  }
  return _closestTriangle != NULL;
}

} // namespace query
} // namespace precice
