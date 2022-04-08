#pragma once

#include <boost/container/flat_map.hpp>

#include "mesh/Mesh.hpp"
#include "precice/types.hpp"

namespace precice {
namespace mesh {

/// Creates a vertex on the destination mesh and adds it to the map to use it to build other connectivity elements
void addVertexToMesh(const Vertex &vertex, boost::container::flat_map<VertexID, Vertex *> &vertexMap, Mesh &destination);

/// Creates an edge on the destination mesh. If withConnection option is true, it also creates the missing vertex for an edge
void addEdgeToMesh(const Edge &edge, boost::container::flat_map<EdgeID, Edge *> &edgeMap, boost::container::flat_map<VertexID, Vertex *> &vertexMap, Mesh &destination, bool withConnection);

/// Creates an face on the destination mesh. If withConnection option is true, it also creates the missing vertices for an edge
void addTriangleToMesh(const Triangle &triangle, boost::container::flat_map<EdgeID, Edge *> &edgeMap, boost::container::flat_map<VertexID, Vertex *> &vertexMap, Mesh &destination, bool withConnection);

/** filters the source Mesh and adds it to the destination Mesh
 * @param[inout] destination the destination mesh to append the filtered Mesh to
 * @param[in] source the source Mesh to filter
 * @param[in] p the filter as a UnaryPredicate on mesh::Vertex 
 */
template <typename UnaryPredicate>
void filterMesh(Mesh &destination, const Mesh &source, UnaryPredicate p, bool withConnection = false)
{
  // Create a flat_map which can contain all vertices of the original mesh.
  // This prevents resizes during the map build-up.
  boost::container::flat_map<VertexID, Vertex *> vertexMap;
  vertexMap.reserve(source.vertices().size());

  for (const Vertex &vertex : source.vertices()) {
    if (p(vertex)) {
      addVertexToMesh(vertex, vertexMap, destination);
    }
  }

  // Create a flat_map which can contain all edges of the original mesh.
  // This prevents resizes during the map build-up.
  boost::container::flat_map<EdgeID, Edge *> edgeMap;
  edgeMap.reserve(source.edges().size());

  // Add all edges formed by the contributing vertices
  for (const Edge &edge : source.edges()) {
    addEdgeToMesh(edge, edgeMap, vertexMap, destination, withConnection);
  }

  // Add all triangles formed by the contributing edges
  if (source.getDimensions() == 3) {
    for (const Triangle &triangle : source.triangles()) {
      addTriangleToMesh(triangle, edgeMap, vertexMap, destination, withConnection);
    }
  }
}

} // namespace mesh
} // namespace precice
