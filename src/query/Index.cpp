#include <Eigen/Core>
#include <algorithm>
#include <boost/iterator/function_output_iterator.hpp>
#include <boost/range/irange.hpp>
#include <utility>

#include "logging/LogMacros.hpp"
#include "precice/types.hpp"
#include "query/Index.hpp"
#include "query/impl/RTreeAdapter.hpp"
#include "utils/Event.hpp"

namespace precice {
extern bool syncMode;
namespace query {

precice::logging::Logger Index::_log{"query::Index"};

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

using VertexTraits   = impl::RTreeTraits<mesh::Vertex>;
using EdgeTraits     = impl::RTreeTraits<mesh::Edge>;
using TriangleTraits = impl::RTreeTraits<mesh::Triangle>;

struct MeshIndices {
  VertexTraits::Ptr   vertexRTree;
  EdgeTraits::Ptr     edgeRTree;
  TriangleTraits::Ptr triangleRTree;
};

class Index::IndexImpl {
public:
  VertexTraits::Ptr   getVertexRTree(const mesh::Mesh &mesh);
  EdgeTraits::Ptr     getEdgeRTree(const mesh::Mesh &mesh);
  TriangleTraits::Ptr getTriangleRTree(const mesh::Mesh &mesh);

  void clear();

private:
  MeshIndices indices;
};

VertexTraits::Ptr Index::IndexImpl::getVertexRTree(const mesh::Mesh &mesh)
{
  if (indices.vertexRTree) {
    return indices.vertexRTree;
  }

  precice::utils::Event e("query.index.getVertexIndexTree." + mesh.getName());

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  impl::RTreeParameters     params;
  VertexTraits::IndexGetter ind(mesh.vertices());
  auto                      tree = std::make_shared<VertexTraits::RTree>(
      boost::irange<std::size_t>(0lu, mesh.vertices().size()), params, ind);

  indices.vertexRTree = std::move(tree);
  return indices.vertexRTree;
}

EdgeTraits::Ptr Index::IndexImpl::getEdgeRTree(const mesh::Mesh &mesh)
{
  if (indices.edgeRTree) {
    return indices.edgeRTree;
  }

  precice::utils::Event e("query.index.getEdgeIndexTree." + mesh.getName());

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  impl::RTreeParameters   params;
  EdgeTraits::IndexGetter ind(mesh.edges());

  auto tree = std::make_shared<EdgeTraits::RTree>(
      boost::irange<std::size_t>(0lu, mesh.edges().size()), params, ind);

  indices.edgeRTree = std::move(tree);
  return indices.edgeRTree;
}

TriangleTraits::Ptr Index::IndexImpl::getTriangleRTree(const mesh::Mesh &mesh)
{
  if (indices.triangleRTree) {
    return indices.triangleRTree;
  }

  precice::utils::Event e("query.index.getTriangleIndexTree." + mesh.getName());

  // We first generate the values for the triangle rtree.
  // The resulting vector is a random access range, which can be passed to the
  // constructor of the rtree for more efficient indexing.
  std::vector<TriangleTraits::IndexType> elements;
  elements.reserve(mesh.triangles().size());
  for (size_t i = 0; i < mesh.triangles().size(); ++i) {
    auto box = bg::return_envelope<RTreeBox>(mesh.triangles()[i]);
    elements.emplace_back(std::move(box), i);
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do.
  impl::RTreeParameters       params;
  TriangleTraits::IndexGetter ind;

  auto tree             = std::make_shared<TriangleTraits::RTree>(elements, params, ind);
  indices.triangleRTree = std::move(tree);
  return indices.triangleRTree;
}

void Index::IndexImpl::clear()
{
  indices.vertexRTree.reset();
  indices.edgeRTree.reset();
  indices.triangleRTree.reset();
}

//
// query::Index
//

Index::Index(mesh::PtrMesh mesh)
    : _mesh(mesh.get())
{
  _pimpl = std::make_unique<IndexImpl>(IndexImpl{});
}

Index::Index(mesh::Mesh &mesh)
    : _mesh(&mesh)
{
  _pimpl = std::make_unique<IndexImpl>(IndexImpl{});
}

// Required for the pimpl idiom to work with std::unique_ptr
Index::~Index() = default;

VertexMatch Index::getClosestVertex(const Eigen::VectorXd &sourceCoord)
{
  PRECICE_TRACE();

  PRECICE_ASSERT(not _mesh->vertices().empty(), _mesh->getName());
  VertexMatch match;
  const auto &rtree = _pimpl->getVertexRTree(*_mesh);
  rtree->query(bgi::nearest(sourceCoord, 1), boost::make_function_output_iterator([&](size_t matchID) {
                 match = VertexMatch(bg::distance(sourceCoord, _mesh->vertices()[matchID]), matchID);
               }));
  return match;
}

std::vector<EdgeMatch> Index::getClosestEdges(const Eigen::VectorXd &sourceCoord, int n)
{
  PRECICE_TRACE();

  const auto &rtree = _pimpl->getEdgeRTree(*_mesh);

  std::vector<EdgeMatch> matches;
  rtree->query(bgi::nearest(sourceCoord, n), boost::make_function_output_iterator([&](size_t matchID) {
                 matches.emplace_back(bg::distance(sourceCoord, _mesh->edges()[matchID]), matchID);
               }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<TriangleMatch> Index::getClosestTriangles(const Eigen::VectorXd &sourceCoord, int n)
{
  PRECICE_TRACE();
  const auto &rtree = _pimpl->getTriangleRTree(*_mesh);

  std::vector<TriangleMatch> matches;
  rtree->query(bgi::nearest(sourceCoord, n),
               boost::make_function_output_iterator([&](TriangleTraits::IndexType const &match) {
                 matches.emplace_back(bg::distance(sourceCoord, _mesh->triangles()[match.second]), match.second);
               }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<VertexID> Index::getVerticesInsideBox(const mesh::Vertex &centerVertex, double radius)
{
  PRECICE_TRACE();

  // Prepare boost::geometry box
  auto coords    = centerVertex.getCoords();
  auto searchBox = query::makeBox(coords.array() - radius, coords.array() + radius);

  const auto &          rtree = _pimpl->getVertexRTree(*_mesh);
  std::vector<VertexID> matches;
  rtree->query(bgi::intersects(searchBox) and bg::index::satisfies([&](size_t const i) { return bg::distance(centerVertex, _mesh->vertices()[i]) <= radius; }),
               std::back_inserter(matches));
  return matches;
}

std::vector<VertexID> Index::getVerticesInsideBox(const mesh::BoundingBox &bb)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  const auto &          rtree = _pimpl->getVertexRTree(*_mesh);
  std::vector<VertexID> matches;
  rtree->query(bgi::intersects(query::makeBox(bb.minCorner(), bb.maxCorner())), std::back_inserter(matches));
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

ProjectionMatch Index::findCellInterpolation(const Eigen::VectorXd &location, int n)
{
  if (_mesh->getDimensions() == 2) {
    auto matchedTriangles = getClosestTriangles(location, n);
    for (const auto &match : matchedTriangles) {
      auto polation = mapping::Polation(location, _mesh->triangles()[match.index]);
      if (polation.isInterpolation()) {
        return {polation, 0.0};
      }
    }

    // If no triangle is found, fall-back on NN
    return findVertexProjection(location);
  } else {
    PRECICE_ASSERT(_mesh->getDimensions() == 3, "Volume coupling 3D not  implemented");
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

void Index::clear()
{
  _pimpl->clear();
}

} // namespace query
} // namespace precice
