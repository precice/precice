#pragma once

#include <vector>
#include "mesh/Edge.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace query {

/// Struct to hold index and distance information of the closest primitive
template <class Tag>
struct MatchType {
  double distance;
  int    index;

  MatchType() = default;
  MatchType(double d, int i)
      : distance(d), index(i){};
  constexpr bool operator<(MatchType const &other) const
  {
    return distance < other.distance;
  };
};

/// Match tags for each primitive type
using VertexMatch   = MatchType<struct VertexMatchTag>;
using EdgeMatch     = MatchType<struct EdgeMatchTag>;
using TriangleMatch = MatchType<struct TriangleTag>;

namespace impl {
struct MeshIndices;
}

class Index {

public:
  Index(const mesh::PtrMesh &mesh);

  std::vector<VertexMatch>   getClosestVertex(const mesh::Vertex &sourceVertex, int n = 1);
  std::vector<EdgeMatch>     getClosestEdge(const mesh::Vertex &sourcesVertex, int n = 1);
  std::vector<TriangleMatch> getClosestTriangle(const mesh::Vertex &sourceVertex, int n = 1);
  std::vector<size_t>        getVerticesInsideBox(const mesh::Vertex &centerVertex, double radius);

  static void clearCache();
  static void clearCache(int meshID);

private:
  std::unique_ptr<impl::MeshIndices> _cache;
  const mesh::PtrMesh                _mesh;
};

} // namespace query
} // namespace precice