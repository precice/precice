#include "query/RTree.hpp"
#include <algorithm>
#include <boost/range/irange.hpp>
#include <cstddef>
#include <utility>
#include <vector>
#include "math/barycenter.hpp"
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace query {

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

// Initialize static member
std::map<int, rtree::MeshIndices> precice::query::rtree::_cachedTrees;

rtree::MeshIndices &rtree::cacheEntry(int meshID)
{
  auto result = _cachedTrees.emplace(std::make_pair(meshID, rtree::MeshIndices{}));
  return result.first->second;
}

VertexTraits::Ptr rtree::getVertexRTree(const mesh::PtrMesh &mesh)
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

EdgeTraits::Ptr rtree::getEdgeRTree(const mesh::PtrMesh &mesh)
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

TriangleTraits::Ptr rtree::getTriangleRTree(const mesh::PtrMesh &mesh)
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

std::vector<MatchType> rtree::getClosestVertex(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n)
{
  auto                   tree = getVertexRTree(targetMesh);
  std::vector<MatchType> matches;
  tree->query(bgi::nearest(source, n), boost::make_function_output_iterator([&](size_t matchID) {
                matches.emplace_back(boost::geometry::distance(source, targetMesh->vertices()[matchID]), matchID);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<MatchType> rtree::getClosestEdge(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n)
{
  auto                   tree = getEdgeRTree(targetMesh);
  std::vector<MatchType> matches;
  tree->query(bgi::nearest(source, n), boost::make_function_output_iterator([&](size_t matchID) {
                matches.emplace_back(boost::geometry::distance(source, targetMesh->edges()[matchID]), matchID);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<MatchType> rtree::getClosestTriangle(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n)
{
  auto                   tree = getTriangleRTree(targetMesh);
  std::vector<MatchType> matches;
  tree->query(bgi::nearest(source, n),
              boost::make_function_output_iterator([&](TriangleTraits::IndexType const &match) {
                matches.emplace_back(bg::distance(source, targetMesh->triangles()[match.second]), match.second);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<size_t> rtree::getVerticesInsideBox(const Box3d &searchBox, const mesh::PtrMesh &targetMesh, const mesh::Vertex &centerVertex, double supportRadius)
{
  auto                tree = getVertexRTree(targetMesh);
  std::vector<size_t> matches;
  tree->query(bgi::intersects(searchBox) and bg::index::satisfies([&](size_t const i) { return bg::distance(centerVertex, targetMesh->vertices()[i]) <= supportRadius; }),
              std::back_inserter(matches));
  return matches;
}

void rtree::tagAllInsideBox(const mesh::BoundingBox &boundingBox, const mesh::PtrMesh &targetMesh)
{
  auto tree = getVertexRTree(targetMesh);
  tree->query(bgi::intersects(RTreeBox(boundingBox.minCorner(), boundingBox.maxCorner())),
              boost::make_function_output_iterator([&targetMesh](size_t idx) {
                targetMesh->vertices()[idx].tag();
              }));
}

void rtree::clear(mesh::Mesh &mesh)
{
  _cachedTrees.erase(mesh.getID());
}

void rtree::clear()
{
  _cachedTrees.clear();
}

Box3d getEnclosingBox(mesh::Vertex const &middlePoint, double sphereRadius)
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
