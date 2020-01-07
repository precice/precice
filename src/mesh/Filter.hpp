#pragma once

#include <boost/container/flat_map.hpp>
#include "mesh/Mesh.hpp"

namespace precice {
namespace mesh {

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
  boost::container::flat_map<int, Vertex *> vertexMap;
  vertexMap.reserve(source.vertices().size());

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

  // Create a flat_map which can contain all edges of the original mesh.
  // This prevents resizes during the map build-up.
  boost::container::flat_map<int, Edge *> edgeMap;
  edgeMap.reserve(source.edges().size());

  // Add all edges formed by the contributing vertices
  for (const Edge &edge : source.edges()) {
    int vertexIndex1 = edge.vertex(0).getID();
    int vertexIndex2 = edge.vertex(1).getID();
    if (vertexMap.count(vertexIndex1) == 1 &&
        vertexMap.count(vertexIndex2) == 1) {
      Edge &e               = destination.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
      edgeMap[edge.getID()] = &e;
    }
  }

  // Add all triangles formed by the contributing edges
  if (source.getDimensions() == 3) {
    for (const Triangle &triangle : source.triangles()) {
      int edgeIndex1 = triangle.edge(0).getID();
      int edgeIndex2 = triangle.edge(1).getID();
      int edgeIndex3 = triangle.edge(2).getID();
      if (edgeMap.count(edgeIndex1) == 1 &&
          edgeMap.count(edgeIndex2) == 1 &&
          edgeMap.count(edgeIndex3) == 1) {
        destination.createTriangle(*edgeMap[edgeIndex1], *edgeMap[edgeIndex2], *edgeMap[edgeIndex3]);
      }
    }
  }
}

} // namespace mesh
} // namespace precice
