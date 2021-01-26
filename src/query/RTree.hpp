#pragma once

#include <vector>
#include "mesh/BoundingBox.hpp"
#include "mesh/Edge.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace query {

class MeshIndices;

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

class rtree {
public:
  /// Get the closest vertex/edge/triangle to a vertex, return indices and distance
  static std::vector<MatchType> getClosestVertex(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n = 1);
  static std::vector<MatchType> getClosestEdge(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n = 1);
  static std::vector<MatchType> getClosestTriangle(const mesh::Vertex &source, const mesh::PtrMesh &targetMesh, int n = 1);

  /// Get all the vertices inside the box around a vertex and radius and surrounded by support radius
  static std::vector<size_t> getVerticesInsideBox(const mesh::PtrMesh &targetMesh, const mesh::Vertex &centerVertex, double supportRadius);

  /// Tag all the vertices which are inside the box which is around the vertex and the radius
  static void tagAllInsideBox(const mesh::BoundingBox &searchBox, const mesh::PtrMesh &targetMesh);

  /// Only clear the trees of that specific mesh
  static void clear(mesh::Mesh &mesh);

  /// Clear the complete cache
  static void clear();
};

} // namespace query

} // namespace precice
