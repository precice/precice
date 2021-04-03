#pragma once

#include "mesh/Filter.hpp"

namespace precice {
namespace mesh {

void addVertexToMesh(const Vertex &vertex, boost::container::flat_map<int, Vertex *> &vertexMap, Mesh &destination)
{
  Vertex &v = destination.createVertex(vertex.getCoords());
  v.setGlobalIndex(vertex.getGlobalIndex());
  if (vertex.isTagged()) {
    v.tag();
  }
  v.setOwner(vertex.isOwner());
  vertexMap[vertex.getID()] = &v;
}

void addEdgeToMesh(const Edge &edge, boost::container::flat_map<int, Edge *> &edgeMap, boost::container::flat_map<int, Vertex *> &vertexMap, Mesh &destination)
{
  if (vertexMap.count(edge.vertex(0).getID()) == 0) {
    addVertexToMesh(edge.vertex(0), vertexMap, destination);
  }
  if (vertexMap.count(edge.vertex(1).getID()) == 0) {
    addVertexToMesh(edge.vertex(1), vertexMap, destination);
  }
  Edge &e               = destination.createEdge(*vertexMap[edge.vertex(0).getID()], *vertexMap[edge.vertex(1).getID()]);
  edgeMap[edge.getID()] = &e;
}

void addTriangleToMesh(const Triangle &triangle, boost::container::flat_map<int, Edge *> &edgeMap, boost::container::flat_map<int, Vertex *> &vertexMap, Mesh &destination)
{
  if (edgeMap.count(triangle.edge(0).getID()) == 0) {
    addEdgeToMesh(triangle.edge(0), edgeMap, vertexMap, destination);
  }
  if (edgeMap.count(triangle.edge(1).getID()) == 0) {
    addEdgeToMesh(triangle.edge(1), edgeMap, vertexMap, destination);
  }
  if (edgeMap.count(triangle.edge(2).getID()) == 0) {
    addEdgeToMesh(triangle.edge(2), edgeMap, vertexMap, destination);
  }
  Triangle &t = destination.createTriangle(*edgeMap[triangle.edge(0).getID()], *edgeMap[triangle.edge(1).getID()], *edgeMap[triangle.edge(2).getID()]);
}

} // namespace mesh
} // namespace precice
