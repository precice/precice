#pragma once

#include <memory>
#include <vector>

#include "logging/Logger.hpp"
#include "mapping/Polation.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Edge.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Tetrahedron.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "precice/types.hpp"

namespace precice {
namespace query {

/// Type used for the IDs of matching entities
using MatchID = int;

constexpr MatchID NO_MATCH{-1};

/// Type used for the distance of a match to the queried point
using Distance = double;
constexpr double INVALID_DISTANCE{-1};

/// Struct to hold the index of a primitive match
template <class Tag>
struct MatchType {
  MatchID index{NO_MATCH};
  MatchType() = default;
  explicit MatchType(MatchID i)
      : index(i){};
};

/// Match tags for each primitive type
using GenericMatch  = MatchType<struct GenericMatchTag>;
using VertexMatch   = MatchType<struct VertexMatchTag>;
using EdgeMatch     = MatchType<struct EdgeMatchTag>;
using TriangleMatch = MatchType<struct TriangleTag>;
using TetraMatch    = MatchType<struct TetraTag>;

/// Struct representing a projection match
struct ProjectionMatch {
  ProjectionMatch(const mapping::Polation &p)
      : polation(p) {}
  ProjectionMatch(mapping::Polation &&p)
      : polation(std::move(p)) {}

  ProjectionMatch(const ProjectionMatch &other) = default;
  ProjectionMatch(ProjectionMatch &&other)      = default;
  ProjectionMatch &operator=(const ProjectionMatch &other) = default;
  ProjectionMatch &operator=(ProjectionMatch &&other) = default;

  mapping::Polation polation;

  bool operator<(ProjectionMatch const &other) const
  {
    return polation.distance() < other.polation.distance();
  };
};

/// Class to query the index trees of the mesh
class Index {

public:
  Index(mesh::PtrMesh mesh);
  Index(mesh::Mesh &mesh);
  ~Index();

  /// Get the closest vertex to the given vertex
  VertexMatch getClosestVertex(const Eigen::VectorXd &sourceCoord);

  /// Get n number of closest vertices to the given vertex
  std::vector<VertexID> getClosestVertices(const Eigen::VectorXd &sourceCoord, int n);

  /// Get n number of closest edges to the given vertex
  std::vector<EdgeMatch> getClosestEdges(const Eigen::VectorXd &sourceCoord, int n);

  /// Get n number of closest triangles to the given vertex
  std::vector<TriangleMatch> getClosestTriangles(const Eigen::VectorXd &sourceCoord, int n);

  /// Return all the vertices inside the box formed by vertex and radius (boundary exclusive)
  std::vector<VertexID> getVerticesInsideBox(const mesh::Vertex &centerVertex, double radius);

  /// Return all the vertices inside a bounding box
  std::vector<VertexID> getVerticesInsideBox(const mesh::BoundingBox &bb);

  /// Return all the tetrahedra whose axis-aligned bounding box contains a vertex
  std::vector<TetrahedronID> getEnclosingTetrahedra(const Eigen::VectorXd &location);

  /**
   * @brief Find the closest interpolation element to the given location.
   * If exists, triangle or edge projection element is returned. If not vertex projection element, which is the nearest neighbor is returned.
   *
   * param[in] sourceVertex
   * param[in] n how many nearest edges/faces are going to be checked
   *
   * param[out] pair of interpolation and the distance to corresponding vertex/edge/triangle
   *
   */
  ProjectionMatch findNearestProjection(const Eigen::VectorXd &location, int n);

  ProjectionMatch findCellOrProjection(const Eigen::VectorXd &location, int n);

  // Index tree, bounds
  mesh::BoundingBox getRtreeBounds();

  /// Clear the index
  void clear();

private:
  class IndexImpl;
  std::unique_ptr<IndexImpl> _pimpl;

  /// The indexed Mesh.
  mesh::Mesh *_mesh;

  static precice::logging::Logger _log;

  /// Closest vertex projection element is always the nearest neighbor
  ProjectionMatch findVertexProjection(const Eigen::VectorXd &location);

  /// Find closest edge interpolation element. If cannot be found, it falls back to vertex projection
  ProjectionMatch findEdgeProjection(const Eigen::VectorXd &location, int n, ProjectionMatch closestVertex);

  /// Find closest face interpolation element. If cannot be found, it falls back to first edge interpolation element, then vertex if necessary
  ProjectionMatch findTriangleProjection(const Eigen::VectorXd &location, int n, ProjectionMatch closestVertex);
};

} // namespace query
} // namespace precice
