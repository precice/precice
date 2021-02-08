#include "Index.hpp"
#include <boost/range/irange.hpp>
#include "impl/Indexer.hpp"
#include "logging/LogMacros.hpp"
#include "utils/Event.hpp"

namespace precice {
extern bool syncMode;
namespace query {

precice::logging::Logger Index::_log{"query::Index"};

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

Index::Index(const mesh::PtrMesh &mesh)
    : _mesh(mesh)
{
  _cache = std::make_unique<impl::MeshIndices>(impl::MeshIndices{});
}

Index::~Index()
{
}

VertexMatch Index::getClosestVertex(const mesh::Vertex &sourceVertex)
{
  PRECICE_TRACE();
  precice::utils::Event event("query.index.getClosestVerticesOnMesh." + _mesh->getName(), precice::syncMode);
  // Add tree to the local cache
  if (not _cache->vertexRTree) {
    _cache->vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
  }

  std::vector<VertexMatch> matches;
  _cache->vertexRTree->query(bgi::nearest(sourceVertex, 1), boost::make_function_output_iterator([&](size_t matchID) {
                               matches.emplace_back(bg::distance(sourceVertex, _mesh->vertices()[matchID]), matchID);
                             }));
  event.stop();
  return matches.back();
}

std::vector<VertexMatch> Index::getClosestVertices(const mesh::Vertex &sourceVertex, int n)
{
  PRECICE_TRACE();
  precice::utils::Event event("query.index.getClosestVerticesOnMesh." + _mesh->getName(), precice::syncMode);
  // Add tree to the local cache
  if (not _cache->vertexRTree) {
    _cache->vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
  }

  std::vector<VertexMatch> matches;
  _cache->vertexRTree->query(bgi::nearest(sourceVertex, n), boost::make_function_output_iterator([&](size_t matchID) {
                               matches.emplace_back(bg::distance(sourceVertex, _mesh->vertices()[matchID]), matchID);
                             }));
  std::sort(matches.begin(), matches.end());
  event.stop();
  return matches;
}

std::vector<EdgeMatch> Index::getClosestEdges(const mesh::Vertex &sourceVertex, int n)
{
  PRECICE_TRACE();
  precice::utils::Event event("query.index.getClosestEdgesOnMesh." + _mesh->getName(), precice::syncMode);
  // Add tree to the local cache
  if (not _cache->edgeRTree) {
    _cache->edgeRTree = impl::Indexer::instance()->getEdgeRTree(_mesh);
  }

  std::vector<EdgeMatch> matches;
  _cache->edgeRTree->query(bgi::nearest(sourceVertex, n), boost::make_function_output_iterator([&](size_t matchID) {
                             matches.emplace_back(bg::distance(sourceVertex, _mesh->edges()[matchID]), matchID);
                           }));
  std::sort(matches.begin(), matches.end());
  event.stop();
  return matches;
}

std::vector<TriangleMatch> Index::getClosestTriangles(const mesh::Vertex &sourceVertex, int n)
{
  PRECICE_TRACE();
  precice::utils::Event event("query.index.getClosestTrianglesOnMesh." + _mesh->getName(), precice::syncMode);
  // Add tree to the local cache
  if (not _cache->triangleRTree) {
    _cache->triangleRTree = impl::Indexer::instance()->getTriangleRTree(_mesh);
  }

  std::vector<TriangleMatch> matches;
  _cache->triangleRTree->query(bgi::nearest(sourceVertex, n),
                               boost::make_function_output_iterator([&](impl::TriangleTraits::IndexType const &match) {
                                 matches.emplace_back(bg::distance(sourceVertex, _mesh->triangles()[match.second]), match.second);
                               }));
  std::sort(matches.begin(), matches.end());
  event.stop();
  return matches;
}

std::vector<size_t> Index::getVerticesInsideBox(const mesh::Vertex &centerVertex, double radius)
{
  PRECICE_TRACE();
  precice::utils::Event event("query.index.getVerticesInsideBoxOnMesh." + _mesh->getName(), precice::syncMode);
  // Add tree to the local cache
  if (_cache->vertexRTree == nullptr) {
    _cache->vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
  }

  auto                searchBox = impl::Indexer::instance()->getEnclosingBox(centerVertex, radius);
  std::vector<size_t> matches;
  _cache->vertexRTree->query(bgi::intersects(searchBox) and bg::index::satisfies([&](size_t const i) { return bg::distance(centerVertex, _mesh->vertices()[i]) <= radius; }),
                             std::back_inserter(matches));
  event.stop();
  return matches;
}

std::vector<size_t> Index::getVerticesInsideBox(const mesh::BoundingBox &bb)
{
  PRECICE_TRACE();
  precice::utils::Event event("query.index.getVerticesInsideBoxOnMesh." + _mesh->getName(), precice::syncMode);
  // Add tree to the local cache
  if (not _cache->vertexRTree) {
    _cache->vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
  }
  std::vector<size_t> matches;
  _cache->vertexRTree->query(bgi::intersects(query::RTreeBox{bb.minCorner(), bb.maxCorner()}), std::back_inserter(matches));
  event.stop();
  return matches;
}

void rtree::clearCache()
{
  impl::Indexer::instance()->clearCache();
}

void rtree::clearCache(int meshID)
{
  impl::Indexer::instance()->clearCache(meshID);
}

void rtree::clearCache(mesh::Mesh &mesh)
{
  impl::Indexer::instance()->clearCache(mesh.getID());
}

} // namespace query
} // namespace precice