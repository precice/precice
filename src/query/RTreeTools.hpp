#pragma once

#include <map>
#include "impl/RTreeAdapter.hpp"
#include "mesh/BoundingBox.hpp"

namespace precice {
namespace query {

/// Class to seperate boost::geometry
class RTreeTools {
public:
  /// Returns either already cached vertex RTree or creates the RTree
  static VertexTraits::Ptr getVertexRTree(const mesh::PtrMesh &mesh);

  /// Returns either already cached edge RTree or creates the RTree
  static EdgeTraits::Ptr getEdgeRTree(const mesh::PtrMesh &mesh);

  /// Returns either already cached triangle RTree or creates the RTree
  static TriangleTraits::Ptr getTriangleRTree(const mesh::PtrMesh &mesh);

  /// Returns the RTree cache for the given mesh
  static MeshIndices &cacheEntry(int MeshID);

  /// Returns a boost::geometry box that encloses a sphere of given radius around a middle point
  static Box3d getEnclosingBox(mesh::Vertex const &middlePoint, double sphereRadius);

  /// Only clear the trees of that specific mesh
  static void clear(mesh::Mesh &mesh);

  /// Clear the complete cache
  static void clear();

  friend struct testing::accessors::rtree;

private:
  /// Cache for all index trees
  static std::map<int, query::MeshIndices> _cachedTrees;
};

} // namespace query

namespace testing {
namespace accessors {
struct rtree {
  static std::map<int, query::MeshIndices> &getCache()
  {
    return query::RTreeTools::_cachedTrees;
  }
};
} // namespace accessors
} // namespace testing

} // namespace precice