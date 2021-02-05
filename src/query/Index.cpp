#include "Index.hpp"
#include <boost/range/irange.hpp>
#include "impl/RTreeWrapper.hpp"
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

std::vector<VertexMatch> Index::getClosestVertex(const mesh::Vertex &sourceVertex, int n)
{
  PRECICE_TRACE();
  precice::utils::Event event("query.index.getClosestVerticesOnMesh." + _mesh->getName(), precice::syncMode);
  // Add tree to the local cache
  if (not _cache->vertexRTree) {
    _cache->vertexRTree = impl::RTreeWrapper::getVertexRTree(_mesh);
  }

  std::vector<VertexMatch> matches;
  _cache->vertexRTree->query(bgi::nearest(sourceVertex, n), boost::make_function_output_iterator([&](size_t matchID) {
                               matches.emplace_back(bg::distance(sourceVertex, _mesh->vertices()[matchID]), matchID);
                             }));
  std::sort(matches.begin(), matches.end());
  event.stop();
  return matches;
}

std::vector<EdgeMatch> Index::getClosestEdge(const mesh::Vertex &sourceVertex, int n)
{
  PRECICE_TRACE();
  precice::utils::Event event("query.index.getClosestEdgesOnMesh." + _mesh->getName(), precice::syncMode);
  // Add tree to the local cache
  if (not _cache->edgeRTree) {
    _cache->edgeRTree = impl::RTreeWrapper::getEdgeRTree(_mesh);
  }

  std::vector<EdgeMatch> matches;
  _cache->edgeRTree->query(bgi::nearest(sourceVertex, n), boost::make_function_output_iterator([&](size_t matchID) {
                             matches.emplace_back(bg::distance(sourceVertex, _mesh->edges()[matchID]), matchID);
                           }));
  std::sort(matches.begin(), matches.end());
  event.stop();
  return matches;
}

std::vector<TriangleMatch> Index::getClosestTriangle(const mesh::Vertex &sourceVertex, int n)
{
  PRECICE_TRACE();
  precice::utils::Event event("query.index.getClosestTrianglesOnMesh." + _mesh->getName(), precice::syncMode);
  // Add tree to the local cache
  if (not _cache->triangleRTree) {
    _cache->triangleRTree = impl::RTreeWrapper::getTriangleRTree(_mesh);
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
    _cache->vertexRTree = impl::RTreeWrapper::getVertexRTree(_mesh);
  }

  auto                searchBox = impl::RTreeWrapper::getEnclosingBox(centerVertex, radius);
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
    _cache->vertexRTree = impl::RTreeWrapper::getVertexRTree(_mesh);
  }
  std::vector<size_t> matches;
  _cache->vertexRTree->query(bgi::intersects(query::RTreeBox{bb.minCorner(), bb.maxCorner()}), std::back_inserter(matches));
  event.stop();
  return matches;
}

void Index::clearCache()
{
  impl::RTreeWrapper::clearCache();
}

void Index::clearCache(int meshID)
{
  impl::RTreeWrapper::clearCache(meshID);
}

void Index::clearCache(mesh::Mesh &mesh)
{
  impl::RTreeWrapper::clearCache(mesh.getID());
}

} // namespace query
} // namespace precice