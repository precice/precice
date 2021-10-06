#include <boost/range/irange.hpp>

#include "Indexer.hpp"
#include "mesh/BoundingBox.hpp"
#include "precice/types.hpp"

namespace precice {
namespace query {
namespace impl {

std::shared_ptr<Indexer> Indexer::instance()
{
  static std::shared_ptr<Indexer> indexer{new Indexer};
  return indexer;
}

MeshIndices &Indexer::cacheEntry(MeshID meshID)
{
  auto result = _cachedTrees.emplace(std::make_pair(meshID, MeshIndices{}));
  return result.first->second;
}

VertexTraits::Ptr Indexer::getVertexRTree(const mesh::PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.vertexRTree) {
    return cache.vertexRTree;
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  impl::RTreeParameters     params;
  VertexTraits::IndexGetter ind(mesh->vertices());
  auto                      tree = std::make_shared<VertexTraits::RTree>(
      boost::irange<std::size_t>(0lu, mesh->vertices().size()), params, ind);

  cache.vertexRTree = std::move(tree);
  return cache.vertexRTree;
}

EdgeTraits::Ptr Indexer::getEdgeRTree(const mesh::PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.edgeRTree) {
    return cache.edgeRTree;
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  impl::RTreeParameters   params;
  EdgeTraits::IndexGetter ind(mesh->edges());
  auto                    tree = std::make_shared<EdgeTraits::RTree>(
      boost::irange<std::size_t>(0lu, mesh->edges().size()), params, ind);

  cache.edgeRTree = std::move(tree);
  return cache.edgeRTree;
}

TriangleTraits::Ptr Indexer::getTriangleRTree(const mesh::PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.triangleRTree) {
    return cache.triangleRTree;
  }

  // We first generate the values for the triangle rtree.
  // The resulting vector is a random access range, which can be passed to the
  // constructor of the rtree for more efficient indexing.
  std::vector<TriangleTraits::IndexType> elements;
  elements.reserve(mesh->triangles().size());
  for (size_t i = 0; i < mesh->triangles().size(); ++i) {
    auto box = bg::return_envelope<RTreeBox>(mesh->triangles()[i]);
    elements.emplace_back(std::move(box), i);
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do.
  impl::RTreeParameters       params;
  TriangleTraits::IndexGetter ind;
  auto                        tree = std::make_shared<TriangleTraits::RTree>(elements, params, ind);
  cache.triangleRTree              = std::move(tree);
  return cache.triangleRTree;
}

size_t Indexer::getCacheSize()
{
  return _cachedTrees.size();
}

void Indexer::clearCache()
{
  _cachedTrees.clear();
}

void Indexer::clearCache(MeshID meshID)
{
  _cachedTrees.erase(meshID);
}

} // namespace impl
} // namespace query
} // namespace precice
