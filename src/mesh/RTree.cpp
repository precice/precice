#include "mesh/RTree.hpp"
#include <algorithm>
#include <boost/range/irange.hpp>
#include <cstddef>
#include <utility>
#include <vector>
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {

namespace bg = boost::geometry;

// Initialize static member
std::map<int, rtree::MeshIndices> precice::mesh::rtree::_cached_trees;

rtree::MeshIndices &rtree::cacheEntry(int meshID)
{
  auto result = _cached_trees.emplace(std::make_pair(meshID, rtree::MeshIndices{}));
  return result.first->second;
}

rtree::vertex_traits::Ptr rtree::getVertexRTree(const PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.vertices) {
    return cache.vertices;
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  RTreeParameters            params;
  vertex_traits::IndexGetter ind(mesh->vertices());
  auto                       tree = std::make_shared<vertex_traits::RTree>(
      boost::irange<std::size_t>(0lu, mesh->vertices().size()), params, ind);

  cache.vertices = tree;
  return tree;
}

rtree::edge_traits::Ptr rtree::getEdgeRTree(const PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.edges) {
    return cache.edges;
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  RTreeParameters          params;
  edge_traits::IndexGetter ind(mesh->edges());
  auto                     tree = std::make_shared<edge_traits::RTree>(
      boost::irange<std::size_t>(0lu, mesh->edges().size()), params, ind);

  cache.edges = tree;
  return tree;
}

rtree::triangle_traits::Ptr rtree::getTriangleRTree(const PtrMesh &mesh)
{
  PRECICE_ASSERT(mesh);
  auto &cache = cacheEntry(mesh->getID());
  if (cache.triangles) {
    return cache.triangles;
  }

  // We first generate the values for the triangle rtree.
  // The resulting vector is a random access range, which can be passed to the
  // constructor of the rtree for more efficient indexing.
  std::vector<triangle_traits::IndexType> elements;
  elements.reserve(mesh->triangles().size());
  for (size_t i = 0; i < mesh->triangles().size(); ++i) {
    auto box = bg::return_envelope<RTreeBox>(mesh->triangles()[i]);
    elements.emplace_back(std::move(box), i);
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do.
  RTreeParameters              params;
  triangle_traits::IndexGetter ind;
  auto                         tree = std::make_shared<triangle_traits::RTree>(elements, params, ind);
  cache.triangles                   = tree;
  return tree;
}

void rtree::clear(Mesh &mesh)
{
  _cached_trees.erase(mesh.getID());
}

void rtree::clear()
{
  _cached_trees.clear();
}

Box3d getEnclosingBox(Vertex const &middlePoint, double sphereRadius)
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

} // namespace mesh
} // namespace precice
