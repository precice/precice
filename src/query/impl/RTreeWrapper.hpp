#pragma once

#include <map>
#include "RTreeAdapter.hpp"

namespace precice {

namespace testing {
namespace accessors {
struct rtree;
}
} // namespace testing

namespace query {
namespace impl {

using VertexTraits   = impl::RTreeTraits<mesh::Vertex>;
using EdgeTraits     = impl::RTreeTraits<mesh::Edge>;
using TriangleTraits = impl::RTreeTraits<mesh::Triangle>;

struct MeshIndices {
  VertexTraits::Ptr   vertexRTree;
  EdgeTraits::Ptr     edgeRTree;
  TriangleTraits::Ptr triangleRTree;
};

/// Class to encapsulate boost::geometry implementations
class RTreeWrapper {
public:
  /// Return vertex index tree from cache, if cache is empty, create the tree
  VertexTraits::Ptr getVertexRTree(const mesh::PtrMesh &mesh);

  /// Return edge index tree from cache, if cache is empty, create the tree
  EdgeTraits::Ptr getEdgeRTree(const mesh::PtrMesh &mesh);

  /// Return triangle index tree from cache, if cache is empty, create the tree
  TriangleTraits::Ptr getTriangleRTree(const mesh::PtrMesh &mesh);

  /// Return boost::geometry version of enclosing box around a point
  Box3d getEnclosingBox(const mesh::Vertex &middlePoint, double sphereRadius);

  /// Clear the whole cache
  static void clearCache();

  /// Clear the cache only for the given mesh
  static void clearCache(int meshID);

private:
  friend struct testing::accessors::rtree;

  MeshIndices &                     cacheEntry(int meshID);
  static std::map<int, MeshIndices> _cachedTrees;
};

} // namespace impl
} // namespace query

namespace testing {
namespace accessors {
struct rtree {
  static std::map<int, query::impl::MeshIndices> &getCache()
  {
    return query::impl::RTreeWrapper::_cachedTrees;
  }
};

} // namespace accessors
} // namespace testing

} // namespace precice