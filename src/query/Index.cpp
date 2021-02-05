#include "Index.hpp"
#include <boost/range/irange.hpp>
#include "impl/RTreeWrapper.hpp"

namespace precice {
namespace query {

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

Index::Index(const mesh::PtrMesh &mesh)
    : _mesh(mesh)
{
}

Index::~Index()
{
}

std::vector<VertexMatch> Index::getClosestVertex(const mesh::Vertex &sourceVertex, int n)
{
  auto tree = impl::RTreeWrapper::getVertexRTree(_mesh);

  std::vector<VertexMatch> matches;
  tree->query(bgi::nearest(sourceVertex, n), boost::make_function_output_iterator([&](size_t matchID) {
                matches.emplace_back(bg::distance(sourceVertex, _mesh->vertices()[matchID]), matchID);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<EdgeMatch> Index::getClosestEdge(const mesh::Vertex &sourceVertex, int n)
{
  auto tree = impl::RTreeWrapper::getEdgeRTree(_mesh);

  std::vector<EdgeMatch> matches;
  tree->query(bgi::nearest(sourceVertex, n), boost::make_function_output_iterator([&](size_t matchID) {
                matches.emplace_back(bg::distance(sourceVertex, _mesh->edges()[matchID]), matchID);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<TriangleMatch> Index::getClosestTriangle(const mesh::Vertex &sourceVertex, int n)
{
  auto tree = impl::RTreeWrapper::getTriangleRTree(_mesh);

  std::vector<TriangleMatch> matches;
  tree->query(bgi::nearest(sourceVertex, n),
              boost::make_function_output_iterator([&](impl::TriangleTraits::IndexType const &match) {
                matches.emplace_back(bg::distance(sourceVertex, _mesh->triangles()[match.second]), match.second);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<size_t> Index::getVerticesInsideBox(const mesh::Vertex &centerVertex, double radius)
{
  auto tree = impl::RTreeWrapper::getVertexRTree(_mesh);

  auto                searchBox = impl::RTreeWrapper::getEnclosingBox(centerVertex, radius);
  std::vector<size_t> matches;
  tree->query(bgi::intersects(searchBox) and bg::index::satisfies([&](size_t const i) { return bg::distance(centerVertex, _mesh->vertices()[i]) <= radius; }),
              std::back_inserter(matches));
  return matches;
}

std::vector<size_t> Index::getVerticesInsideBox(const mesh::BoundingBox &bb)
{
  auto                tree = impl::RTreeWrapper::getVertexRTree(_mesh);
  std::vector<size_t> matches;
  tree->query(bgi::intersects(query::RTreeBox{bb.minCorner(), bb.maxCorner()}), std::back_inserter(matches));
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