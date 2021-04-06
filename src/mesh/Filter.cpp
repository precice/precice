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

void addEdgeToMesh(const Edge &edge, boost::container::flat_map<int, Edge *> &edgeMap, boost::container::flat_map<int, Vertex *> &vertexMap, Mesh &destination, bool withConnection)
{
  int vertexIndex1 = edge.vertex(0).getID();
  int vertexIndex2 = edge.vertex(1).getID();

  if (vertexMap.count(vertexIndex1) == 0 and
      vertexMap.count(vertexIndex2) == 0) {
    return;
  }

  if (withConnection) {
    if (vertexMap.count(vertexIndex1) == 0) {
      addVertexToMesh(edge.vertex(0), vertexMap, destination);
    }
    if (vertexMap.count(vertexIndex2) == 0) {
      addVertexToMesh(edge.vertex(1), vertexMap, destination);
    }
  }

  if (vertexMap.count(vertexIndex1) == 1 and
      vertexMap.count(vertexIndex2) == 1) {
    Edge &e               = destination.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
    edgeMap[edge.getID()] = &e;
  }
}

void addTriangleToMesh(const Triangle &triangle, boost::container::flat_map<int, Edge *> &edgeMap, boost::container::flat_map<int, Vertex *> &vertexMap, Mesh &destination, bool withConnection)
{
  int edgeIndex1 = triangle.edge(0).getID();
  int edgeIndex2 = triangle.edge(1).getID();
  int edgeIndex3 = triangle.edge(2).getID();

  if (edgeMap.count(edgeIndex1) == 0 and
      edgeMap.count(edgeIndex2) == 0 and
      edgeMap.count(edgeIndex3) == 0) {
    return;
  }

  if (withConnection) {
    if (edgeMap.count(edgeIndex1) == 0) {
      addEdgeToMesh(triangle.edge(0), edgeMap, vertexMap, destination, withConnection);
    }
    if (edgeMap.count(edgeIndex2) == 0) {
      addEdgeToMesh(triangle.edge(1), edgeMap, vertexMap, destination, withConnection);
    }
    if (edgeMap.count(edgeIndex3) == 0) {
      addEdgeToMesh(triangle.edge(2), edgeMap, vertexMap, destination, withConnection);
    }
  }
  if (edgeMap.count(edgeIndex1) == 1 and
      edgeMap.count(edgeIndex2) == 1 and
      edgeMap.count(edgeIndex3) == 1) {
    Triangle &t = destination.createTriangle(*edgeMap[edgeIndex1], *edgeMap[edgeIndex2], *edgeMap[edgeIndex3]);
  }
}

} // namespace mesh
} // namespace precice
