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

std::vector<VertexMatch> Index::getClosestVertex(const mesh::Vertex &sourceVertex, int n)
{

  const auto &tree = impl::RTreeWrapper::getVertexRTree(_mesh);

  std::vector<VertexMatch> matches;
  tree->query(bgi::nearest(sourceVertex, n), boost::make_function_output_iterator([&](size_t matchID) {
                matches.emplace_back(bg::distance(sourceVertex, _mesh->vertices()[matchID]), matchID);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<EdgeMatch> Index::getClosestEdge(const mesh::Vertex &sourceVertex, int n)
{

  const auto &tree = impl::RTreeWrapper::getEdgeRTree(_mesh);

  std::vector<EdgeMatch> matches;
  tree->query(bgi::nearest(sourceVertex, n), boost::make_function_output_iterator([&](size_t matchID) {
                matches.emplace_back(bg::distance(sourceVertex, _mesh->edges()[matchID]), matchID);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

std::vector<TriangleMatch> Index::getClosestTriangle(const mesh::Vertex &sourceVertex, int n)
{
  const auto &tree = impl::RTreeWrapper::getTriangleRTree(_mesh);

  std::vector<TriangleMatch> matches;
  tree->query(bgi::nearest(sourceVertex, n),
              boost::make_function_output_iterator([&](impl::TriangleTraits::IndexType const &match) {
                matches.emplace_back(bg::distance(sourceVertex, _mesh->triangles()[match.second]), match.second);
              }));
  std::sort(matches.begin(), matches.end());
  return matches;
}

} // namespace query
} // namespace precice