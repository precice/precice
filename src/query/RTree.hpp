#pragma once

#include <Eigen/Core>
#include <boost/geometry.hpp>
#include <boost/range/irange.hpp>
#include <iosfwd>
#include <map>
#include <memory>
#include <type_traits>
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "query/impl/RTreeAdapter.hpp"
#include "utils/Statistics.hpp"

namespace precice {
namespace testing {
namespace accessors {
class rtree;
} // namespace accessors
} // namespace testing
} // namespace precice

namespace precice {
namespace query {

class rtree {
public:
  /// Returns the pointer to boost::geometry::rtree for the given mesh vertices: Right now only used by tests
  /*
   * Creates and fills the tree, if it wasn't requested before, otherwise it returns the cached tree.
   */
  static VertexTraits::Ptr   getVertexRTree(const mesh::PtrMesh &mesh);
  static EdgeTraits::Ptr     getEdgeRTree(const mesh::PtrMesh &mesh);
  static TriangleTraits::Ptr getTriangleRTree(const mesh::PtrMesh &mesh);

  // MAYBE ENUMED CALLS
  /// Get the closest vertex/edge/triangle to a vertex, return indices and distance
  static std::vector<MatchType> getClosestVertex(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n = 1);
  static std::vector<MatchType> getClosestEdge(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n = 1);
  static std::vector<MatchType> getClosestTriangle(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n = 1);

  /// Get all the vertices inside the box and surrounded by support radius
  static std::vector<size_t> getVerticesInsideBox(const Box3d &searchBox, const mesh::PtrMesh &targetMesh, const mesh::Vertex &centerVertex, double supportRadius);

  /// Tag all the vertices which are inside the given box
  static void tagAllInsideBox(const mesh::BoundingBox &searchBox, const mesh::PtrMesh &targetMesh);

  /// Only clear the trees of that specific mesh
  static void clear(mesh::Mesh &mesh);

  /// Clear the complete cache
  static void clear();

  friend struct testing::accessors::rtree;

private:
  // REMOVE THESE FROM INTERFACE/HEADERS
  struct MeshIndices {
    VertexTraits::Ptr   vertices;
    EdgeTraits::Ptr     edges;
    TriangleTraits::Ptr triangles;
  };

  static MeshIndices &cacheEntry(int MeshID);

  using RTreeCache = std::map<int, MeshIndices>;
  static RTreeCache _cachedTrees; ///< Cache for all index trees
};

/// Returns a boost::geometry box that encloses a sphere of given radius around a middle point
Box3d getEnclosingBox(mesh::Vertex const &middlePoint, double sphereRadius);

} // namespace query

namespace testing {
namespace accessors {
struct rtree {
  static query::rtree::RTreeCache &getCache()
  {
    return query::rtree::_cachedTrees;
  }
};
} // namespace accessors
} // namespace testing

} // namespace precice
