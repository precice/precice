#pragma once

#include <map>
#include "RTreeAdapter.hpp"

namespace precice {
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

class RTreeWrapper {
public:
  static VertexTraits::Ptr   getVertexRTree(const mesh::PtrMesh &mesh);
  static EdgeTraits::Ptr     getEdgeRTree(const mesh::PtrMesh &mesh);
  static TriangleTraits::Ptr getTriangleRTree(const mesh::PtrMesh &mesh);

private:
  static MeshIndices &              cacheEntry(int meshID);
  static std::map<int, MeshIndices> _cachedTrees;
};

} // namespace impl
} // namespace query
} // namespace precice