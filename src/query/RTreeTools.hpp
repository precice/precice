#pragma once

#include <map>
#include "impl/RTreeAdapter.hpp"
#include "mesh/BoundingBox.hpp"

namespace precice {
namespace query {

class RTreeTools {
public:
  static VertexTraits::Ptr   getVertexRTree(const mesh::PtrMesh &mesh);
  static EdgeTraits::Ptr     getEdgeRTree(const mesh::PtrMesh &mesh);
  static TriangleTraits::Ptr getTriangleRTree(const mesh::PtrMesh &mesh);

  /// Returns a boost::geometry box that encloses a sphere of given radius around a middle point
  static Box3d getEnclosingBox(mesh::Vertex const &middlePoint, double sphereRadius);

  /// Only clear the trees of that specific mesh
  static void clear(mesh::Mesh &mesh);

  /// Clear the complete cache
  static void clear();

  friend struct testing::accessors::rtree;

private:
  static MeshIndices &                     cacheEntry(int MeshID);
  static std::map<int, query::MeshIndices> _cachedTrees; ///< Cache for all index trees
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