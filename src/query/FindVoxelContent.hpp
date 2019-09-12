#pragma once

#include <Eigen/Core>
#include "mesh/Group.hpp"

namespace precice {
namespace query {

/// Find visitables that are contained in a given voxel.
class FindVoxelContent {
public:
  /// Deciders to include/exclude elements touching the voxels bounds.
  enum BoundaryInclusion {
    INCLUDE_BOUNDARY,
    EXCLUDE_BOUNDARY
  };

  /**
   * @brief Constructor.
   *
   * @param[in] voxelCenter       Coordinates of the voxel center
   * @param[in] halflengths       Half of the sidelengths of the voxel
   * @param[in] boundaryInclusion Decides if boundaries should be included or excluded
   */
  FindVoxelContent(
      const Eigen::VectorXd &voxelCenter,
      const Eigen::VectorXd &halflengths,
      BoundaryInclusion      boundaryInclusion);

  /**
   * @brief Performs the find operation on the content of the container.
   *
   * In two dimensions, the content edges need to have valid normals. In three
   * dimensions, the triangles need to have valid normals.
   */
  template <typename CONTAINER_T>
  bool operator()(CONTAINER_T &container);

  /// Returns voxel center coordinates
  const Eigen::VectorXd &getVoxelCenter() const;

  /// Returns halflength of voxel
  const Eigen::VectorXd &getVoxelHalflengths() const;

  /// Returns true, if touching objects are included in content
  BoundaryInclusion getBoundaryInclusion() const;

  /// Returns a Group with all contained visitables
  mesh::Group &content();

  void clear();

private:
  mutable logging::Logger _log{"query::FindVoxelContent"};

  /// Center point of the voxel
  Eigen::VectorXd _voxelCenter;

  /// Half the sidelengths of the voxel
  Eigen::VectorXd _voxelHalflengths;

  /// Determines, whether objects that touch the voxel are contained
  BoundaryInclusion _boundaryInclusion;

  /// Spatial dimensions of the voxel (2D or 3D).
  int _dimensions;

  /// Contained and partially contained visitables of the voxel
  mesh::Group _content;

  void checkVertex(mesh::Vertex &vertex);

  void checkEdge(mesh::Edge &edge);

  void checkTriangle(mesh::Triangle &triangle);

  /// Returns true, if a plane square and a segment intersect.
  bool computeIntersection(
      const Eigen::Vector3d &squareCenter,
      const Eigen::VectorXd &halflengths,
      int                    squareNormalDirection,
      const Eigen::Vector3d &firstPointSegment,
      const Eigen::Vector3d &secondPointSegment,
      bool                   countTouchingAsIntersection) const;

  bool computeIntersection(
      const mesh::Triangle & triangle,
      const Eigen::Vector3d &firstPointSegment,
      const Eigen::Vector3d &secondPointSegment,
      bool                   countTouchingAsIntersection);
};

// ---------------------------------------------------------- HEAER DEFINITIONS

template <typename CONTAINER_T>
bool FindVoxelContent::operator()(
    CONTAINER_T &container)
{
  for (mesh::Vertex &vertex : container.vertices()) {
    checkVertex(vertex);
  }
  for (mesh::Edge &edge : container.edges()) {
    checkEdge(edge);
  }
  for (mesh::Triangle &triangle : container.triangles()) {
    checkTriangle(triangle);
  }
  return not _content.empty();
}

} // namespace query
} // namespace precice
