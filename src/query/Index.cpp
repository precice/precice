#include "Index.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
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

struct Index::IndexImpl {
  impl::MeshIndices indices;
};

Index::Index(const mesh::PtrMesh &mesh)
    : _mesh(mesh)
{
  _pimpl = std::make_unique<IndexImpl>(IndexImpl{});
}

Index::~Index()
{
}

VertexMatch Index::getClosestVertex(const Eigen::VectorXd &sourceCoord)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  if (not _pimpl->indices.vertexRTree) {
    precice::utils::Event event("query.index.getVertexIndexTree." + _mesh->getName(), precice::syncMode);
    _pimpl->indices.vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
    event.stop();
  }

  VertexMatch match;
  _pimpl->indices.vertexRTree->query(bgi::nearest(sourceCoord, 1), boost::make_function_output_iterator([&](size_t matchID) {
                                       match = VertexMatch(bg::distance(sourceCoord, _mesh->vertices()[matchID]), matchID);
                                     }));
  return match;
}

std::vector<EdgeMatch> Index::getClosestEdges(const Eigen::VectorXd &sourceCoord, int n)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  if (not _pimpl->indices.edgeRTree) {
    precice::utils::Event event("query.index.getEdgeIndexTree." + _mesh->getName(), precice::syncMode);
    _pimpl->indices.edgeRTree = impl::Indexer::instance()->getEdgeRTree(_mesh);
    event.stop();
  }

  std::vector<EdgeMatch> matches;
  _pimpl->indices.edgeRTree->query(bgi::nearest(sourceCoord, n), boost::make_function_output_iterator([&](size_t matchID) {
                                     matches.emplace_back(bg::distance(sourceCoord, _mesh->edges()[matchID]), matchID);
                                   }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<TriangleMatch> Index::getClosestTriangles(const Eigen::VectorXd &sourceCoord, int n)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  if (not _pimpl->indices.triangleRTree) {
    precice::utils::Event event("query.index.getTriangleIndexTree." + _mesh->getName(), precice::syncMode);
    _pimpl->indices.triangleRTree = impl::Indexer::instance()->getTriangleRTree(_mesh);
    event.stop();
  }

  std::vector<TriangleMatch> matches;
  _pimpl->indices.triangleRTree->query(bgi::nearest(sourceCoord, n),
                                       boost::make_function_output_iterator([&](impl::TriangleTraits::IndexType const &match) {
                                         matches.emplace_back(bg::distance(sourceCoord, _mesh->triangles()[match.second]), match.second);
                                       }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<size_t> Index::getVerticesInsideBox(const mesh::Vertex &centerVertex, double radius)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  if (_pimpl->indices.vertexRTree == nullptr) {
    precice::utils::Event event("query.index.getVertexIndexTree." + _mesh->getName(), precice::syncMode);
    _pimpl->indices.vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
    event.stop();
  }

  // Prepare boost::geometry box
  auto coords    = centerVertex.getCoords();
  auto searchBox = query::makeBox(coords.array() - radius, coords.array() + radius);

  std::vector<size_t> matches;
  _pimpl->indices.vertexRTree->query(bgi::intersects(searchBox) and bg::index::satisfies([&](size_t const i) { return bg::distance(centerVertex, _mesh->vertices()[i]) <= radius; }),
                                     std::back_inserter(matches));
  return matches;
}

std::vector<size_t> Index::getVerticesInsideBox(const mesh::BoundingBox &bb)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  if (not _pimpl->indices.vertexRTree) {
    precice::utils::Event event("query.index.getVertexIndexTree." + _mesh->getName(), precice::syncMode);
    _pimpl->indices.vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
    event.stop();
  }
  std::vector<size_t> matches;
  _pimpl->indices.vertexRTree->query(bgi::intersects(query::makeBox(bb.minCorner(), bb.maxCorner())), std::back_inserter(matches));
  return matches;
}

std::pair<mapping::Polation, double> Index::findNearestProjection(const Eigen::VectorXd &location, int n)
{
  if (_mesh->getDimensions() == 2) {
    return findEdgeProjection(location, n);
  } else {
    return findTriangleProjection(location, n);
  }
}

std::pair<mapping::Polation, double> Index::findVertexProjection(const Eigen::VectorXd &location)
{
  auto match = getClosestVertex(location);
  return std::pair<mapping::Polation, double>(mapping::Polation(_mesh->vertices()[match.index]), match.distance);
}

std::pair<mapping::Polation, double> Index::findEdgeProjection(const Eigen::VectorXd &location, int n)
{
  auto matchedEdges = getClosestEdges(location, n);
  for (const auto &match : matchedEdges) {
    auto polation = mapping::Polation(location, _mesh->edges()[match.index]);
    if (polation.isInterpolation()) {
      return std::pair<mapping::Polation, double>(polation, match.distance);
    }
  }
  // Could not find edge projection element, fall back to vertex projection
  return findVertexProjection(location);
}

std::pair<mapping::Polation, double> Index::findTriangleProjection(const Eigen::VectorXd &location, int n)
{
  auto matchedTriangles = getClosestTriangles(location, n);
  for (const auto &match : matchedTriangles) {
    auto polation = mapping::Polation(location, _mesh->triangles()[match.index]);
    if (polation.isInterpolation()) {
      return std::pair<mapping::Polation, double>(polation, match.distance);
    }
  }

  // Could not triangle find projection element, fall back to edge projection
  return findEdgeProjection(location, n);
}

void clearCache()
{
  impl::Indexer::instance()->clearCache();
}

void clearCache(int meshID)
{
  impl::Indexer::instance()->clearCache(meshID);
}

void clearCache(mesh::Mesh &mesh)
{
  impl::Indexer::instance()->clearCache(mesh.getID());
}

} // namespace query
} // namespace precice
