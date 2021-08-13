#include <Eigen/Core>
#include <algorithm>
#include <boost/range/irange.hpp>
#include <utility>

#include "Index.hpp"
#include "impl/Indexer.hpp"
#include "logging/LogMacros.hpp"
#include "precice/types.hpp"
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

Index::Index(mesh::PtrMesh mesh)
    : _mesh(std::move(mesh))
{
  _pimpl = std::make_unique<IndexImpl>(IndexImpl{});
}

Index::~Index() = default;

VertexMatch Index::getClosestVertex(const Eigen::VectorXd &sourceCoord)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  if (not _pimpl->indices.vertexRTree) {
    precice::utils::Event e("query.index.getVertexIndexTree." + _mesh->getName());
    _pimpl->indices.vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
  }

  PRECICE_ASSERT(not _mesh->vertices().empty(), _mesh->getName());
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
    precice::utils::Event e("query.index.getEdgeIndexTree." + _mesh->getName());
    _pimpl->indices.edgeRTree = impl::Indexer::instance()->getEdgeRTree(_mesh);
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
    precice::utils::Event e("query.index.getTriangleIndexTree." + _mesh->getName());
    _pimpl->indices.triangleRTree = impl::Indexer::instance()->getTriangleRTree(_mesh);
  }

  std::vector<TriangleMatch> matches;
  _pimpl->indices.triangleRTree->query(bgi::nearest(sourceCoord, n),
                                       boost::make_function_output_iterator([&](impl::TriangleTraits::IndexType const &match) {
                                         matches.emplace_back(bg::distance(sourceCoord, _mesh->triangles()[match.second]), match.second);
                                       }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<VertexID> Index::getVerticesInsideBox(const mesh::Vertex &centerVertex, double radius)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  if (_pimpl->indices.vertexRTree == nullptr) {
    precice::utils::Event e("query.index.getVertexIndexTree." + _mesh->getName());
    _pimpl->indices.vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
  }

  // Prepare boost::geometry box
  auto coords    = centerVertex.getCoords();
  auto searchBox = query::makeBox(coords.array() - radius, coords.array() + radius);

  std::vector<VertexID> matches;
  _pimpl->indices.vertexRTree->query(bgi::intersects(searchBox) and bg::index::satisfies([&](size_t const i) { return bg::distance(centerVertex, _mesh->vertices()[i]) <= radius; }),
                                     std::back_inserter(matches));
  return matches;
}

std::vector<VertexID> Index::getVerticesInsideBox(const mesh::BoundingBox &bb)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  if (not _pimpl->indices.vertexRTree) {
    precice::utils::Event e("query.index.getVertexIndexTree." + _mesh->getName());
    _pimpl->indices.vertexRTree = impl::Indexer::instance()->getVertexRTree(_mesh);
  }
  std::vector<VertexID> matches;
  _pimpl->indices.vertexRTree->query(bgi::intersects(query::makeBox(bb.minCorner(), bb.maxCorner())), std::back_inserter(matches));
  return matches;
}

ProjectionMatch Index::findNearestProjection(const Eigen::VectorXd &location, int n)
{
  if (_mesh->getDimensions() == 2) {
    return findEdgeProjection(location, n);
  } else {
    return findTriangleProjection(location, n);
  }
}

ProjectionMatch Index::findVertexProjection(const Eigen::VectorXd &location)
{
  auto match = getClosestVertex(location);
  return {mapping::Polation{_mesh->vertices()[match.index]}, match.distance};
}

ProjectionMatch Index::findEdgeProjection(const Eigen::VectorXd &location, int n)
{
  auto matchedEdges = getClosestEdges(location, n);
  for (const auto &match : matchedEdges) {
    auto polation = mapping::Polation(location, _mesh->edges()[match.index]);
    if (polation.isInterpolation()) {
      return {polation, match.distance};
    }
  }
  // Could not find edge projection element, fall back to vertex projection
  return findVertexProjection(location);
}

ProjectionMatch Index::findTriangleProjection(const Eigen::VectorXd &location, int n)
{
  auto matchedTriangles = getClosestTriangles(location, n);
  for (const auto &match : matchedTriangles) {
    auto polation = mapping::Polation(location, _mesh->triangles()[match.index]);
    if (polation.isInterpolation()) {
      return {polation, match.distance};
    }
  }

  // Could not triangle find projection element, fall back to edge projection
  return findEdgeProjection(location, n);
}

void clearCache()
{
  impl::Indexer::instance()->clearCache();
}

void clearCache(MeshID meshID)
{
  impl::Indexer::instance()->clearCache(meshID);
}

void clearCache(mesh::Mesh &mesh)
{
  impl::Indexer::instance()->clearCache(mesh.getID());
}

} // namespace query
} // namespace precice
