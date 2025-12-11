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
void filterMesh(Mesh &destination, const Mesh &source, UnaryPredicate p)
{
  // Create a flat_map which can contain all vertices of the original mesh.
  // This prevents resizes during the map build-up.
  boost::container::flat_map<VertexID, Vertex *> vertexMap;
  vertexMap.reserve(source.nVertices());

  for (const Vertex &vertex : source.vertices()) {
    if (p(vertex)) {
      Vertex &v = destination.createVertex(vertex.getCoords());
      v.setGlobalIndex(vertex.getGlobalIndex());
      if (vertex.isTagged())
        v.tag();
      v.setOwner(vertex.isOwner());
      vertexMap[vertex.getID()] = &v;
    }
  }

  // Add all edges formed by the contributing vertices
  for (const Edge &edge : source.edges()) {
    VertexID vertexIndex1 = edge.vertex(0).getID();
    VertexID vertexIndex2 = edge.vertex(1).getID();
    if (vertexMap.count(vertexIndex1) == 1 &&
        vertexMap.count(vertexIndex2) == 1) {
      destination.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
    }
  }

  // Add all triangles formed by the contributing vertices
  for (const Triangle &triangle : source.triangles()) {
    VertexID vertexIndex1 = triangle.vertex(0).getID();
    VertexID vertexIndex2 = triangle.vertex(1).getID();
    VertexID vertexIndex3 = triangle.vertex(2).getID();
    if (vertexMap.count(vertexIndex1) == 1 &&
        vertexMap.count(vertexIndex2) == 1 &&
        vertexMap.count(vertexIndex3) == 1) {
      destination.createTriangle(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2], *vertexMap[vertexIndex3]);
    }
  }

  // Add all tetrahedra formed by the contributing vertices
  for (const Tetrahedron &tetra : source.tetrahedra()) {
    VertexID vertexIndex1 = tetra.vertex(0).getID();
    VertexID vertexIndex2 = tetra.vertex(1).getID();
    VertexID vertexIndex3 = tetra.vertex(2).getID();
    VertexID vertexIndex4 = tetra.vertex(3).getID();
    if (vertexMap.count(vertexIndex1) == 1 &&
        vertexMap.count(vertexIndex2) == 1 &&
        vertexMap.count(vertexIndex3) == 1 &&
        vertexMap.count(vertexIndex4) == 1) {
      destination.createTetrahedron(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2], *vertexMap[vertexIndex3], *vertexMap[vertexIndex4]);
    }
  }
}

} // namespace precice::mesh
