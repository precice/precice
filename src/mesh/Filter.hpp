#pragma once

#include <boost/container/flat_map.hpp>

#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"

namespace precice::mesh {

/** filters the source Mesh and adds it to the destination Mesh
 * @param[inout] destination the destination mesh to append the filtered Mesh to
 * @param[in] source the source Mesh to filter
 * @param[in] p the filter as a UnaryPredicate on mesh::Vertex
 */
template <typename UnaryPredicate>
void filterMesh(Mesh &destination, const Mesh &source, UnaryPredicate &&p)
{
  if (source.hasConnectivity()) {
    filterMeshWithConnectivity(destination, source, std::forward<UnaryPredicate>(p));
  } else {
    filterMeshWithoutConnectivity(destination, source, std::forward<UnaryPredicate>(p));
  }
}

template <typename UnaryPredicate>
void filterMeshWithConnectivity(Mesh &destination, const Mesh &source, UnaryPredicate p)
{
  // Create a lookup table which can contain all vertices of the original mesh.
  std::vector<Vertex *> vertexMap(source.nVertices(), nullptr);

  for (const Vertex &vertex : source.vertices()) {
    if (p(vertex)) {
      Vertex &v = destination.createVertex(vertex.getCoords());
      v.setGlobalIndex(vertex.getGlobalIndex());
      v.setTagged(vertex.isTagged());
      v.setOwner(vertex.isOwner());
      vertexMap[vertex.getID()] = &v;
    }
  }

  auto fetch = [&vertexMap](int vid) -> Vertex * {
#ifndef NDEBUG
    return vertexMap.at(vid);
#else
    return vertexMap[vid];
#endif
  };

  // Add all edges formed by the contributing vertices
  for (const Edge &edge : source.edges()) {
    auto vertex1 = fetch(edge.vertex(0).getID());
    auto vertex2 = fetch(edge.vertex(1).getID());
    if (vertex1 &&
        vertex2) {
      destination.createEdge(*vertex1, *vertex2);
    }
  }

  // Add all triangles formed by the contributing vertices
  for (const Triangle &triangle : source.triangles()) {
    auto vertex1 = fetch(triangle.vertex(0).getID());
    auto vertex2 = fetch(triangle.vertex(1).getID());
    auto vertex3 = fetch(triangle.vertex(2).getID());
    if (vertex1 &&
        vertex2 &&
        vertex3) {
      destination.createTriangle(*vertex1, *vertex2, *vertex3);
    }
  }

  // Add all tetrahedra formed by the contributing vertices
  for (const Tetrahedron &tetra : source.tetrahedra()) {
    auto vertex1 = fetch(tetra.vertex(0).getID());
    auto vertex2 = fetch(tetra.vertex(1).getID());
    auto vertex3 = fetch(tetra.vertex(2).getID());
    auto vertex4 = fetch(tetra.vertex(3).getID());
    if (vertex1 &&
        vertex2 &&
        vertex3 &&
        vertex4) {
      destination.createTetrahedron(*vertex1, *vertex2, *vertex3, *vertex4);
    }
  }
}

template <typename UnaryPredicate>
void filterMeshWithoutConnectivity(Mesh &destination, const Mesh &source, UnaryPredicate p)
{
  for (const Vertex &vertex : source.vertices()) {
    if (p(vertex)) {
      Vertex &v = destination.createVertex(vertex.getCoords());
      v.setGlobalIndex(vertex.getGlobalIndex());
      v.setTagged(vertex.isTagged());
      v.setOwner(vertex.isOwner());
    }
  }
}

} // namespace precice::mesh
