#pragma once

#include <memory>
#include <vector>
#include "FindClosest.hpp"
#include "logging/Logger.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Edge.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace query {

/// Struct to hold index and distance information of the closest primitive
template <class Tag>
struct MatchType {
  double distance;
  int    index;
  MatchType() = default;
  MatchType(double d, int i)
      : distance(d), index(i){};
  constexpr bool operator<(MatchType const &other) const
  {
    return distance < other.distance;
  };
};

/// Match tags for each primitive type
using GenericMatch  = MatchType<struct GenericMatchTag>;
using VertexMatch   = MatchType<struct VertexMatchTag>;
using EdgeMatch     = MatchType<struct EdgeMatchTag>;
using TriangleMatch = MatchType<struct TriangleTag>;

/// Class to query the index trees of the mesh
class Index {

public:
  Index(const mesh::PtrMesh &mesh);
  ~Index();

  /// Get n number of closest vertices to the given vertex
  VertexMatch getClosestVertex(const mesh::Vertex &sourceVertex);
  VertexMatch getClosestVertex(const Eigen::VectorXd &sourceCoord);

  /// Get n number of closest edges to the given vertex
  std::vector<EdgeMatch> getClosestEdges(const mesh::Vertex &sourceVertex, int n);
  std::vector<EdgeMatch> getClosestEdges(const Eigen::VectorXd &sourceCoord, int n);

  /// Get n number of closest triangles to the given vertex
  std::vector<TriangleMatch> getClosestTriangles(const mesh::Vertex &sourceVertex, int n);
  std::vector<TriangleMatch> getClosestTriangles(const Eigen::VectorXd &sourceCoord, int n);

  /// Return all the vertices inside the box formed by vertex and radius
  std::vector<size_t> getVerticesInsideBox(const mesh::Vertex &centerVertex, double radius);

  /// Return all the vertices inside a bounding box
  std::vector<size_t> getVerticesInsideBox(const mesh::BoundingBox &bb);

  /**
   * @brief Find the closest interpolation element to a vertex. 
   * If exists, triangle or edge projection element is returned. If not vertex projection element, which is the neares neighbor is returned.
   * 
   * param[in] sourceVertex 
   * param[in] n how many nearest edges/faces are going to be checked
   * 
   * param[out] pair of interpolation elements and the distance to corresponding vertex/edge/triangle
   *
  */
  std::pair<InterpolationElements, double> findNearestProjection(const mesh::Vertex &sourceVertex, int n);

private:
  struct IndexImpl;
  std::unique_ptr<IndexImpl> _pimpl;

  const mesh::PtrMesh             _mesh;
  static precice::logging::Logger _log;

  /// Closest vertex projection element is always the nearest neighbor
  std::pair<InterpolationElements, double> findVertexProjection(const mesh::Vertex &sourceVertex);

  /// Find closest edge interpolation element. If cannot be found, it falls back to vertex projection
  std::pair<InterpolationElements, double> findEdgeProjection(const mesh::Vertex &sourceVertex, int n);

  /// Find closest face interpolation element. If cannot be found, it falls back to first edge interpolation element, then vertex if necessary
  std::pair<InterpolationElements, double> findTriangleProjection(const mesh::Vertex &sourceVertex, int n);
};

/// Clear all the cache
void clearCache();

/// Clear the cache of given mesh
void clearCache(int meshID);

/// Clear the cache of given mesh
void clearCache(mesh::Mesh &mesh);

} // namespace query
} // namespace precice