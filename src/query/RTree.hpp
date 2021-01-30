#pragma once

#include <vector>
#include "mesh/BoundingBox.hpp"
#include "mesh/Edge.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace query {

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

using VertexMatch   = MatchType<struct VertexMatchTag>;
using EdgeMatch     = MatchType<struct EdgeMatchTag>;
using TriangleMatch = MatchType<struct TriangleTag>;

void addVertexToRTree(const mesh::Vertex &vertex, int meshID);
void addEdgeToRTree(const mesh::Edge &vertex, int meshID);
void addTriangleToRTree(const mesh::Triangle &vertex, int meshID);

/// Get the closest vertex/edge/triangle to a vertex, return indices and distance
std::vector<VertexMatch>   getClosestVertex(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n = 1);
std::vector<EdgeMatch>     getClosestEdge(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n = 1);
std::vector<TriangleMatch> getClosestTriangle(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n = 1);

/// Get all the vertices inside the box around a vertex and radius and surrounded by support radius
std::vector<size_t> getVerticesInsideBox(const mesh::PtrMesh &targetMesh, const mesh::Vertex &centerVertex, double supportRadius);

/// Tag all the vertices which are inside the box which is around the vertex and the radius
void tagAllInsideBox(const mesh::BoundingBox &searchBox, const mesh::PtrMesh &targetMesh);

/// Only clear the trees of that specific mesh
void clearRTreeCache(mesh::Mesh &mesh);

/// Clear the complete cache
void clearRTreeCache();

} // namespace query
} // namespace precice
