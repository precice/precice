#include "query/RTree.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <boost/geometry.hpp>
#include <boost/range/irange.hpp>
#include <cstddef>
#include <utility>
#include <vector>
#include "math/barycenter.hpp"
#include "mesh/Vertex.hpp"
#include "query/RTreeTools.hpp"
#include "query/impl/RTreeAdapter.hpp"
#include "utils/assertion.hpp"

#include <boost/version.hpp>
#if BOOST_VERSION < 106600
#include <boost/function_output_iterator.hpp>
#else
#include <boost/iterator/function_output_iterator.hpp>
#endif

namespace precice {
namespace query {

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

void addVertexToRTree(const mesh::Vertex &vertex, int meshID)
{
  auto &cache = RTreeTools::cacheEntry(meshID);
  if (not cache.vertices) {
    return;
  } else {
    cache.vertices->insert(vertex.getID());
  }
}

void addEdgeToRTree(const mesh::Edge &edge, int meshID)
{
  auto &cache = RTreeTools::cacheEntry(meshID);
  if (not cache.edges) {
    return;
  } else {
    cache.edges->insert(edge.getID());
  }
}

void addTriangleToRTree(const mesh::Triangle &triangle, int meshID)
{
  auto &cache = RTreeTools::cacheEntry(meshID);
  if (not cache.triangles) {
    return;
  } else {
    cache.triangles->insert(TriangleTraits::IndexType(bg::return_envelope<RTreeBox>(triangle), triangle.getID()));
  }
}

std::vector<VertexMatch> getClosestVertex(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n)
{
  auto                     tree = RTreeTools::getVertexRTree(targetMesh);
  std::vector<VertexMatch> matches;
  tree->query(bgi::nearest(source, n), boost::make_function_output_iterator([&](size_t matchID) {
                matches.emplace_back(boost::geometry::distance(source, targetMesh->vertices()[matchID]), matchID);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<EdgeMatch> getClosestEdge(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n)
{
  auto                   tree = RTreeTools::getEdgeRTree(targetMesh);
  std::vector<EdgeMatch> matches;
  tree->query(bgi::nearest(source, n), boost::make_function_output_iterator([&](size_t matchID) {
                matches.emplace_back(boost::geometry::distance(source, targetMesh->edges()[matchID]), matchID);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<TriangleMatch> getClosestTriangle(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n)
{
  auto                       tree = RTreeTools::getTriangleRTree(targetMesh);
  std::vector<TriangleMatch> matches;
  tree->query(bgi::nearest(source, n),
              boost::make_function_output_iterator([&](TriangleTraits::IndexType const &match) {
                matches.emplace_back(bg::distance(source, targetMesh->triangles()[match.second]), match.second);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<size_t> getVerticesInsideBox(const mesh::PtrMesh &targetMesh, const mesh::Vertex &centerVertex, double supportRadius)
{
  auto                searchBox = RTreeTools::getEnclosingBox(centerVertex, supportRadius);
  auto                tree      = RTreeTools::getVertexRTree(targetMesh);
  std::vector<size_t> matches;
  tree->query(bgi::intersects(searchBox) and bg::index::satisfies([&](size_t const i) { return bg::distance(centerVertex, targetMesh->vertices()[i]) <= supportRadius; }),
              std::back_inserter(matches));
  return matches;
}

void tagAllInsideBox(const mesh::BoundingBox &boundingBox, const mesh::PtrMesh &targetMesh)
{
  auto tree = RTreeTools::getVertexRTree(targetMesh);
  tree->query(bgi::intersects(RTreeBox(boundingBox.minCorner(), boundingBox.maxCorner())),
              boost::make_function_output_iterator([&targetMesh](size_t idx) {
                targetMesh->vertices()[idx].tag();
              }));
}

void clearRTreeCache(mesh::Mesh &mesh)
{
  RTreeTools::clear(mesh);
}

void clearRTreeCache()
{
  RTreeTools::clear();
}

} // namespace query
} // namespace precice
