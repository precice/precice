#pragma once

#include <Eigen/Core>
#include <limits>
#include <stddef.h>
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/// Finds the closest Vertex object hold by Mesh objects.
class FindClosestVertex {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] searchPoint Coordinates of origin of search for closest point.
   */
  explicit FindClosestVertex(const Eigen::VectorXd &searchPoint);

  /**
   * @brief Searches among all Vertex objects hold by the given Mesh object.
   *
   * When called for different meshes, the closest Vertex object in among all
   * the meshes is found.
   */
  template <typename CONTAINER_T>
  bool operator()(CONTAINER_T &container);

  /// Returns the coordinates of the search point.
  const Eigen::VectorXd &getSearchPoint() const;

  bool hasFound() const;

  /**
   * @brief Returns the euclidian distance to the closest vertex found.
   *
   * @pre find() has been called and returned true.
   */
  double getEuclidianDistance();

  /**
   * @brief Returns the closest vertex found.
   *
   * @pre find() has been called and returned true.
   */
  mesh::Vertex &getClosestVertex();

private:
  /// Origin of search for closest vertex.
  Eigen::VectorXd _searchPoint;

  /// Distance to closest Vertex object found.
  double _shortestDistance = std::numeric_limits<double>::max();

  /// Pointer to closest Vertex object found.
  mesh::Vertex *_closestVertex = nullptr;
};

// --------------------------------------------------------- HEADER DEFINITIONS

template <typename CONTAINER_T>
bool FindClosestVertex::operator()(
    CONTAINER_T &container)
{
  Eigen::VectorXd vectorDistance = Eigen::VectorXd::Zero(_searchPoint.size());
  for (mesh::Vertex &vertex : container.vertices()) {
    PRECICE_ASSERT(vertex.getDimensions() == _searchPoint.size(),
                   vertex.getDimensions(), _searchPoint.size());
    vectorDistance = vertex.getCoords();
    vectorDistance -= _searchPoint;
    double distance = vectorDistance.norm();
    if (distance < _shortestDistance) {
      _shortestDistance = distance;
      _closestVertex    = &vertex;
    }
  }
  return _closestVertex != NULL;
}

} // namespace query
} // namespace precice
