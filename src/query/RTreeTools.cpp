#include "RTreeTools.hpp"
#include <boost/geometry.hpp>
#include <boost/range/irange.hpp>

namespace precice {
namespace query {

MeshIndices &RTreeTools::cacheEntry(int meshID)
{
  auto result = _cachedTrees.emplace(std::make_pair(meshID, MeshIndices{}));
  return result.first->second;
}

VertexTraits::Ptr RTreeTools::getVertexRTree(const mesh::PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.vertices) {
    return cache.vertices;
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  RTreeParameters           params;
  VertexTraits::IndexGetter ind(mesh->vertices());
  auto                      tree = std::make_shared<VertexTraits::RTree>(
      boost::irange<std::size_t>(0lu, mesh->vertices().size()), params, ind);

  cache.vertices = tree;
  return tree;
}

EdgeTraits::Ptr RTreeTools::getEdgeRTree(const mesh::PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.edges) {
    return cache.edges;
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  RTreeParameters         params;
  EdgeTraits::IndexGetter ind(mesh->edges());
  auto                    tree = std::make_shared<EdgeTraits::RTree>(
      boost::irange<std::size_t>(0lu, mesh->edges().size()), params, ind);

  cache.edges = tree;
  return tree;
}

TriangleTraits::Ptr RTreeTools::getTriangleRTree(const mesh::PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.triangles) {
    return cache.triangles;
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
  RTreeParameters             params;
  TriangleTraits::IndexGetter ind;
  auto                        tree = std::make_shared<TriangleTraits::RTree>(elements, params, ind);
  cache.triangles                  = tree;
  return tree;
}

void RTreeTools::clear(mesh::Mesh &mesh)
{
  _cachedTrees.erase(mesh.getID());
}

void RTreeTools::clear()
{
  _cachedTrees.clear();
}

Box3d RTreeTools::getEnclosingBox(mesh::Vertex const &middlePoint, double sphereRadius)
{
  namespace bg = boost::geometry;
  auto &coords = middlePoint.getCoords();

  Box3d box;
  bg::set<bg::min_corner, 0>(box, bg::get<0>(coords) - sphereRadius);
  bg::set<bg::min_corner, 1>(box, bg::get<1>(coords) - sphereRadius);
  bg::set<bg::min_corner, 2>(box, bg::get<2>(coords) - sphereRadius);

  bg::set<bg::max_corner, 0>(box, bg::get<0>(coords) + sphereRadius);
  bg::set<bg::max_corner, 1>(box, bg::get<1>(coords) + sphereRadius);
  bg::set<bg::max_corner, 2>(box, bg::get<2>(coords) + sphereRadius);

  return box;
}

} // namespace query
} // namespace precice