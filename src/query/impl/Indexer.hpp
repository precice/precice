#pragma once

#include <map>

#include "precice/types.hpp"
#include "query/impl/RTreeAdapter.hpp"

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

/// Class to encapsulate boost::geometry implementations
class Indexer {
public:
  Indexer(const Indexer &) = delete;
  Indexer &operator=(const Indexer &) = delete;

  static std::shared_ptr<Indexer> instance();

  /// Return vertex index tree from cache, if cache is empty, create the tree
  VertexTraits::Ptr getVertexRTree(const mesh::PtrMesh &mesh);

  /// Return edge index tree from cache, if cache is empty, create the tree
  EdgeTraits::Ptr getEdgeRTree(const mesh::PtrMesh &mesh);

  /// Return triangle index tree from cache, if cache is empty, create the tree
  TriangleTraits::Ptr getTriangleRTree(const mesh::PtrMesh &mesh);

  size_t getCacheSize();

  /// Clear the whole cache
  void clearCache();

  /// Clear the cache only for the given mesh
  void clearCache(MeshID meshID);

private:
  Indexer(){};
  MeshIndices &              cacheEntry(MeshID meshID);
  std::map<int, MeshIndices> _cachedTrees;
};

} // namespace impl
} // namespace query
} // namespace precice
