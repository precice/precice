#include <Eigen/Core>
#include <algorithm>
#include <boost/iterator/function_output_iterator.hpp>
#include <boost/range/irange.hpp>
#include <utility>

#include "logging/LogMacros.hpp"
#include "precice/impl/Types.hpp"
#include "profiling/Event.hpp"
#include "query/Index.hpp"
#include "query/impl/RTreeAdapter.hpp"

namespace precice::query {

precice::logging::Logger Index::_log{"query::Index"};

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

using VertexTraits      = impl::RTreeTraits<mesh::Vertex>;
using EdgeTraits        = impl::RTreeTraits<mesh::Edge>;
using TriangleTraits    = impl::RTreeTraits<mesh::Triangle>;
using TetrahedronTraits = impl::RTreeTraits<mesh::Tetrahedron>;

struct MeshIndices {
  VertexTraits::Ptr      vertexRTree;
  EdgeTraits::Ptr        edgeRTree;
  TriangleTraits::Ptr    triangleRTree;
  TetrahedronTraits::Ptr tetraRTree;
};

class Index::IndexImpl {
public:
  VertexTraits::Ptr      getVertexRTree(const mesh::Mesh &mesh);
  EdgeTraits::Ptr        getEdgeRTree(const mesh::Mesh &mesh);
  TriangleTraits::Ptr    getTriangleRTree(const mesh::Mesh &mesh);
  TetrahedronTraits::Ptr getTetraRTree(const mesh::Mesh &mesh);

  void clear();

private:
  MeshIndices indices;
};

VertexTraits::Ptr Index::IndexImpl::getVertexRTree(const mesh::Mesh &mesh)
{
  if (indices.vertexRTree) {
    return indices.vertexRTree;
  }

  precice::profiling::Event e("query.index.getVertexIndexTree." + mesh.getName());

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do. Even passing an index range instead of calling
  // tree->insert repeatedly is about 10x faster.
  impl::RTreeParameters     params;
  VertexTraits::IndexGetter ind(mesh.vertices());
  auto                      tree = std::make_shared<VertexTraits::RTree>(
      boost::irange<std::size_t>(0lu, mesh.nVertices()), params, ind);

  indices.vertexRTree = std::move(tree);
  return indices.vertexRTree;
}

EdgeTraits::Ptr Index::IndexImpl::getEdgeRTree(const mesh::Mesh &mesh)
{
  if (indices.edgeRTree) {
    return indices.edgeRTree;
  }

  precice::profiling::Event e("query.index.getEdgeIndexTree." + mesh.getName());

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

  precice::profiling::Event e("query.index.getTriangleIndexTree." + mesh.getName());

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

TetrahedronTraits::Ptr Index::IndexImpl::getTetraRTree(const mesh::Mesh &mesh)
{
  if (indices.tetraRTree) {
    return indices.tetraRTree;
  }

  precice::profiling::Event e("query.index.getTetraIndexTree." + mesh.getName());

  // We first generate the values for the tetra rtree.
  // The resulting vector is a random access range, which can be passed to the
  // constructor of the rtree for more efficient indexing.
  std::vector<TetrahedronTraits::IndexType> elements;
  elements.reserve(mesh.tetrahedra().size());
  for (size_t i = 0; i < mesh.tetrahedra().size(); ++i) {
    // We use a custom function to compute the AABB, because
    // bg::return_envelope was designed for polygons.
    auto box = makeBox(mesh.tetrahedra()[i]);
    elements.emplace_back(std::move(box), i);
  }

  // Generating the rtree is expensive, so passing everything in the ctor is
  // the best we can do.
  impl::RTreeParameters          params;
  TetrahedronTraits::IndexGetter ind;

  auto tree          = std::make_shared<TetrahedronTraits::RTree>(elements, params, ind);
  indices.tetraRTree = std::move(tree);
  return indices.tetraRTree;
}

void Index::IndexImpl::clear()
{
  indices.vertexRTree.reset();
  indices.edgeRTree.reset();
  indices.triangleRTree.reset();
  indices.tetraRTree.reset();
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

  PRECICE_ASSERT(not _mesh->empty(), _mesh->getName());
  VertexMatch match;
  const auto &rtree = _pimpl->getVertexRTree(*_mesh);
  rtree->query(bgi::nearest(sourceCoord, 1), boost::make_function_output_iterator([&](size_t matchID) {
                 match = VertexMatch(matchID);
               }));
  return match;
}

std::vector<VertexID> Index::getClosestVertices(const Eigen::VectorXd &sourceCoord, int n)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(!(_mesh->empty()), _mesh->getName());
  std::vector<VertexID> matches;
  const auto           &rtree = _pimpl->getVertexRTree(*_mesh);

  rtree->query(bgi::nearest(sourceCoord, n), boost::make_function_output_iterator([&](size_t matchID) {
                 matches.emplace_back(matchID);
               }));
  return matches;
}

std::vector<EdgeMatch> Index::getClosestEdges(const Eigen::VectorXd &sourceCoord, int n)
{
  PRECICE_TRACE();

  const auto &rtree = _pimpl->getEdgeRTree(*_mesh);

  std::vector<EdgeMatch> matches;
  matches.reserve(n);
  rtree->query(bgi::nearest(sourceCoord, n), boost::make_function_output_iterator([&](size_t matchID) {
                 matches.emplace_back(matchID);
               }));
  return matches;
}

std::vector<TriangleMatch> Index::getClosestTriangles(const Eigen::VectorXd &sourceCoord, int n)
{
  PRECICE_TRACE();
  const auto &rtree = _pimpl->getTriangleRTree(*_mesh);

  std::vector<TriangleMatch> matches;
  matches.reserve(n);
  rtree->query(bgi::nearest(sourceCoord, n),
               boost::make_function_output_iterator([&](TriangleTraits::IndexType const &match) {
                 matches.emplace_back(match.second);
               }));
  return matches;
}

std::vector<VertexID> Index::getVerticesInsideBox(const mesh::Vertex &centerVertex, double radius)
{
  PRECICE_TRACE();

  // Prepare boost::geometry box
  auto coords    = centerVertex.getCoords();
  auto searchBox = query::makeBox(coords.array() - radius, coords.array() + radius);

  const auto           &rtree = _pimpl->getVertexRTree(*_mesh);
  std::vector<VertexID> matches;
  rtree->query(bgi::intersects(searchBox) and bg::index::satisfies([&](size_t const i) { return bg::distance(centerVertex, _mesh->vertex(i)) < radius; }),
               std::back_inserter(matches));
  return matches;
}

bool Index::isAnyVertexInsideBox(const mesh::Vertex &centerVertex, double radius)
{
  PRECICE_TRACE();

  // Prepare boost::geometry box
  auto coords    = centerVertex.getCoords();
  auto searchBox = query::makeBox(coords.array() - radius, coords.array() + radius);

  const auto &rtree = _pimpl->getVertexRTree(*_mesh);

  auto queryIter = rtree->qbegin(bgi::intersects(searchBox) and bg::index::satisfies([&](size_t const i) { return bg::distance(centerVertex, _mesh->vertex(i)) < radius; }));
  bool hasMatch  = queryIter != rtree->qend();
  return hasMatch;
}

std::vector<VertexID> Index::getVerticesInsideBox(const mesh::BoundingBox &bb)
{
  PRECICE_TRACE();
  // Add tree to the local cache
  const auto           &rtree = _pimpl->getVertexRTree(*_mesh);
  std::vector<VertexID> matches;
  rtree->query(bgi::intersects(query::makeBox(bb.minCorner(), bb.maxCorner())), std::back_inserter(matches));
  return matches;
}

std::vector<TetrahedronID> Index::getEnclosingTetrahedra(const Eigen::VectorXd &location)
{
  PRECICE_TRACE();
  const auto &rtree = _pimpl->getTetraRTree(*_mesh);

  std::vector<TetrahedronID> matches;
  rtree->query(bgi::covers(location), boost::make_function_output_iterator([&](TetrahedronTraits::IndexType const &match) {
                 matches.emplace_back(match.second);
               }));
  return matches;
}

ProjectionMatch Index::findNearestProjection(const Eigen::VectorXd &location, int n)
{
  if (_mesh->getDimensions() == 2) {
    return findEdgeProjection(location, n, findVertexProjection(location));
  } else {
    return findTriangleProjection(location, n, findVertexProjection(location));
  }
}

ProjectionMatch Index::findCellOrProjection(const Eigen::VectorXd &location, int n)
{
  if (_mesh->getDimensions() == 2) {
    auto matchedTriangles = getClosestTriangles(location, n);
    for (const auto &match : matchedTriangles) {
      auto polation = mapping::Polation(location, _mesh->triangles()[match.index]);
      if (polation.isInterpolation()) {
        return {std::move(polation)};
      }
    }

    // If no triangle is found, fall-back on NP
    return findNearestProjection(location, n);
  } else {

    // Find correct tetra, or fall back to NP
    auto matchedTetra = getEnclosingTetrahedra(location);
    for (const auto &match : matchedTetra) {
      // Matches are raw indices, not (indices, distance) pairs
      auto polation = mapping::Polation(location, _mesh->tetrahedra()[match]);
      if (polation.isInterpolation()) {
        return {std::move(polation)};
      }
    }
    return findNearestProjection(location, n);
  }
}

ProjectionMatch Index::findVertexProjection(const Eigen::VectorXd &location)
{
  auto match = getClosestVertex(location);
  return {mapping::Polation{location, _mesh->vertex(match.index)}};
}

ProjectionMatch Index::findEdgeProjection(const Eigen::VectorXd &location, int n, ProjectionMatch closestVertex)
{
  std::vector<ProjectionMatch> candidates;
  candidates.reserve(n);
  for (const auto &match : getClosestEdges(location, n)) {
    auto polation = mapping::Polation(location, _mesh->edges()[match.index]);
    if (polation.isInterpolation()) {
      candidates.emplace_back(std::move(polation));
    }
  }

  // Could not find edge projection element, fall back to vertex projection
  if (candidates.empty()) {
    return closestVertex;
  }

  // Prefer the closest vertex if it closer than the closest edge
  auto min = std::min_element(candidates.begin(), candidates.end());
  if (min->polation.distance() > closestVertex.polation.distance()) {
    return closestVertex;
  }
  return *min;
}

ProjectionMatch Index::findTriangleProjection(const Eigen::VectorXd &location, int n, ProjectionMatch closestVertex)
{
  std::vector<ProjectionMatch> candidates;
  candidates.reserve(n);
  for (const auto &match : getClosestTriangles(location, n)) {
    auto polation = mapping::Polation(location, _mesh->triangles()[match.index]);
    if (polation.isInterpolation()) {
      candidates.emplace_back(std::move(polation));
    }
  }

  // Could not find triangle projection element, fall back to edge projection
  if (candidates.empty()) {
    return findEdgeProjection(location, n, std::move(closestVertex));
  }

  // Fallback to edge projection if a vertex is closer than the best triangle match
  auto min = std::min_element(candidates.begin(), candidates.end());
  if (min->polation.distance() > closestVertex.polation.distance()) {
    return findEdgeProjection(location, n, std::move(closestVertex));
  }

  return *min;
}

mesh::BoundingBox Index::getRtreeBounds()
{
  PRECICE_TRACE();
  // if the mesh is empty, we will most likely hit an assertion in the bounding box class
  // therefore, we keep the assert here, but might want to return an empty bounding box in case
  // we want to allow calling this function with empty meshes
  PRECICE_ASSERT(_mesh->nVertices() > 0);

  auto            rtreeBox = _pimpl->getVertexRTree(*_mesh)->bounds();
  int             dim      = _mesh->getDimensions();
  Eigen::VectorXd min(dim), max(dim);

  min[0] = rtreeBox.min_corner().get<0>();
  min[1] = rtreeBox.min_corner().get<1>();
  max[0] = rtreeBox.max_corner().get<0>();
  max[1] = rtreeBox.max_corner().get<1>();

  if (dim > 2) {
    min[2] = rtreeBox.min_corner().get<2>();
    max[2] = rtreeBox.max_corner().get<2>();
  }
  return mesh::BoundingBox{min, max};
}

void Index::clear()
{
  _pimpl->clear();
}

} // namespace precice::query
