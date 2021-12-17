#include "mesh/Filter.hpp"
namespace precice {
namespace mesh {

void addVertexToMesh(const Vertex &vertex, boost::container::flat_map<VertexID, Vertex *> &vertexMap, Mesh &destination)
{
  Vertex &v = destination.createVertex(vertex.getCoords());
  v.setGlobalIndex(vertex.getGlobalIndex());
  if (vertex.isTagged()) {
    v.tag();
  }
  v.setOwner(vertex.isOwner());
  vertexMap[vertex.getID()] = &v;
}

void addEdgeToMesh(const Edge &edge, boost::container::flat_map<EdgeID, Edge *> &edgeMap, boost::container::flat_map<VertexID, Vertex *> &vertexMap, Mesh &destination, bool withConnection)
{
  const int  vertexIndex1  = edge.vertex(0).getID();
  const int  vertexIndex2  = edge.vertex(1).getID();
  const bool vertexMarked1 = vertexMap.count(vertexIndex1) == 1;
  const bool vertexMarked2 = vertexMap.count(vertexIndex2) == 1;
  const bool allMarked     = vertexMarked1 && vertexMarked2;
  const bool noneMarked    = !(vertexMarked1 || vertexMarked2);

  if (noneMarked) {
    return;
  }

  if (withConnection) {
    if (!vertexMarked1) {
      addVertexToMesh(edge.vertex(0), vertexMap, destination);
    }
    if (!vertexMarked2) {
      addVertexToMesh(edge.vertex(1), vertexMap, destination);
    }
  }

  if (allMarked || withConnection) {
    Edge &e               = destination.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
    edgeMap[edge.getID()] = &e;
  }
}

void addTriangleToMesh(const Triangle &triangle, boost::container::flat_map<EdgeID, Edge *> &edgeMap, boost::container::flat_map<VertexID, Vertex *> &vertexMap, Mesh &destination, bool withConnection)
{
  const int  edgeIndex1  = triangle.edge(0).getID();
  const int  edgeIndex2  = triangle.edge(1).getID();
  const int  edgeIndex3  = triangle.edge(2).getID();
  const bool edgeMarked1 = edgeMap.count(edgeIndex1) == 1;
  const bool edgeMarked2 = edgeMap.count(edgeIndex2) == 1;
  const bool edgeMarked3 = edgeMap.count(edgeIndex3) == 1;
  const bool allMarked   = edgeMarked1 && edgeMarked2 && edgeMarked3;
  const bool noneMarked  = !(edgeMarked1 || edgeMarked2 || edgeMarked3);

  if (noneMarked) {
    return;
  }

  if (withConnection) {
    if (!edgeMarked1) {
      addEdgeToMesh(triangle.edge(0), edgeMap, vertexMap, destination, withConnection);
    }
    if (!edgeMarked2) {
      addEdgeToMesh(triangle.edge(1), edgeMap, vertexMap, destination, withConnection);
    }
    if (!edgeMarked3) {
      addEdgeToMesh(triangle.edge(2), edgeMap, vertexMap, destination, withConnection);
    }
  }

  if (allMarked || withConnection) {
    destination.createTriangle(*edgeMap[edgeIndex1], *edgeMap[edgeIndex2], *edgeMap[edgeIndex3]);
  }
}

} // namespace mesh
} // namespace precice
