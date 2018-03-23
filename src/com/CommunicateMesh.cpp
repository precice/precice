#include "CommunicateMesh.hpp"
#include <map>
#include <vector>
#include "Communication.hpp"
#include "com/SharedPointer.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"

namespace precice
{
namespace com
{

logging::Logger CommunicateMesh::_log("com::CommunicateMesh");

CommunicateMesh::CommunicateMesh(
    com::PtrCommunication communication)
    : _communication(communication)
{
}

void CommunicateMesh::sendMesh(
    const mesh::Mesh &mesh,
    int               rankReceiver)
{
  TRACE(mesh.getName(), rankReceiver);
  int dim = mesh.getDimensions();

  int numberOfVertices = mesh.vertices().size();
  _communication->send(numberOfVertices, rankReceiver);
  if (numberOfVertices > 0) {
    double coords[numberOfVertices * dim];
    int    globalIDs[numberOfVertices];
    for (int i = 0; i < numberOfVertices; i++) {
      for (int d = 0; d < dim; d++) {
        coords[i * dim + d] = mesh.vertices()[i].getCoords()[d];
      }
      globalIDs[i] = mesh.vertices()[i].getGlobalIndex();
    }
    _communication->send(coords, numberOfVertices * dim, rankReceiver);
    _communication->send(globalIDs, numberOfVertices, rankReceiver);
  }

  int numberOfEdges = mesh.edges().size();
  _communication->send(numberOfEdges, rankReceiver);
  if (numberOfEdges > 0) {
    //we need to send the vertexIDs first such that the right edges can be created later
    //contrary to the normal sendMesh, this variant must also work for adding delta meshes
    int vertexIDs[numberOfVertices];
    for (int i = 0; i < numberOfVertices; i++) {
      vertexIDs[i] = mesh.vertices()[i].getID();
    }
    _communication->send(vertexIDs, numberOfVertices, rankReceiver);

    int edgeIDs[numberOfEdges * 2];
    for (int i = 0; i < numberOfEdges; i++) {
      edgeIDs[i * 2]     = mesh.edges()[i].vertex(0).getID();
      edgeIDs[i * 2 + 1] = mesh.edges()[i].vertex(1).getID();
    }
    _communication->send(edgeIDs, numberOfEdges * 2, rankReceiver);
  }

  if (dim == 3) {
    int numberOfTriangles = mesh.triangles().size();
    _communication->send(numberOfTriangles, rankReceiver);
    if (numberOfTriangles > 0) {
      //we need to send the edgeIDs first such that the right edges can be created later
      //contrary to the normal sendMesh, this variant must also work for adding delta meshes
      int edgeIDs[numberOfEdges];
      for (int i = 0; i < numberOfEdges; i++) {
        edgeIDs[i] = mesh.edges()[i].getID();
      }
      _communication->send(edgeIDs, numberOfEdges, rankReceiver);

      int triangleIDs[numberOfTriangles * 3];
      for (int i = 0; i < numberOfTriangles; i++) {
        triangleIDs[i * 3]     = mesh.triangles()[i].edge(0).getID();
        triangleIDs[i * 3 + 1] = mesh.triangles()[i].edge(1).getID();
        triangleIDs[i * 3 + 2] = mesh.triangles()[i].edge(2).getID();
      }
      _communication->send(triangleIDs, numberOfTriangles * 3, rankReceiver);
    }
  }
}

void CommunicateMesh::receiveMesh(
    mesh::Mesh &mesh,
    int         rankSender)
{
  TRACE(mesh.getName(), rankSender);
  int dim = mesh.getDimensions();

  std::vector<mesh::Vertex *>   vertices;
  std::map<int, mesh::Vertex *> vertexMap;
  int                           numberOfVertices = 0;
  _communication->receive(numberOfVertices, rankSender);
  DEBUG("Number of vertices to receive: " << numberOfVertices);

  if (numberOfVertices > 0) {
    double vertexCoords[numberOfVertices * dim];
    int    globalIDs[numberOfVertices];
    _communication->receive(vertexCoords, numberOfVertices * dim, rankSender);
    _communication->receive(globalIDs, numberOfVertices, rankSender);
    for (int i = 0; i < numberOfVertices; i++) {
      Eigen::VectorXd coords(dim);
      for (int d = 0; d < dim; d++) {
        coords[d] = vertexCoords[i * dim + d];
      }
      mesh::Vertex &v = mesh.createVertex(coords);
      assertion(v.getID() >= 0, v.getID());
      v.setGlobalIndex(globalIDs[i]);
      vertices.push_back(&v);
    }
  }

  int                       numberOfEdges = 0;
  std::vector<mesh::Edge *> edges;
  _communication->receive(numberOfEdges, rankSender);
  DEBUG("Number of edges to receive: " << numberOfEdges);
  if (numberOfEdges > 0) {
    int vertexIDs[numberOfVertices];
    _communication->receive(vertexIDs, numberOfVertices, rankSender);
    for (int i = 0; i < numberOfVertices; i++) {
      vertexMap[vertexIDs[i]] = vertices[i];
    }

    int edgeIDs[numberOfEdges * 2];
    _communication->receive(edgeIDs, numberOfEdges * 2, rankSender);
    for (int i = 0; i < numberOfEdges; i++) {
      assertion(vertexMap.find(edgeIDs[i * 2]) != vertexMap.end());
      assertion(vertexMap.find(edgeIDs[i * 2 + 1]) != vertexMap.end());
      assertion(edgeIDs[i * 2] != edgeIDs[i * 2 + 1]);
      mesh::Edge &e = mesh.createEdge(*vertexMap[edgeIDs[i * 2]], *vertexMap[edgeIDs[i * 2 + 1]]);
      edges.push_back(&e);
    }
  }

  if (dim == 3) {
    int numberOfTriangles = 0;
    _communication->receive(numberOfTriangles, rankSender);
    DEBUG("Number of Triangles to receive: " << numberOfTriangles);
    DEBUG("Number of Edges: " << edges.size());
    if (numberOfTriangles > 0) {
      assertion((edges.size() > 0) || (numberOfTriangles == 0));
      int edgeIDs[numberOfEdges];
      _communication->receive(edgeIDs, numberOfEdges, rankSender);
      std::map<int, mesh::Edge *> edgeMap;
      for (int i = 0; i < numberOfEdges; i++) {
        edgeMap[edgeIDs[i]] = edges[i];
      }

      int triangleIDs[numberOfTriangles * 3];
      _communication->receive(triangleIDs, numberOfTriangles * 3, rankSender);

      for (int i = 0; i < numberOfTriangles; i++) {
        assertion(edgeMap.find(triangleIDs[i * 3]) != edgeMap.end());
        assertion(edgeMap.find(triangleIDs[i * 3 + 1]) != edgeMap.end());
        assertion(edgeMap.find(triangleIDs[i * 3 + 2]) != edgeMap.end());
        assertion(triangleIDs[i * 3] != triangleIDs[i * 3 + 1]);
        assertion(triangleIDs[i * 3 + 1] != triangleIDs[i * 3 + 2]);
        assertion(triangleIDs[i * 3 + 2] != triangleIDs[i * 3]);
        mesh.createTriangle(*edgeMap[triangleIDs[i * 3]], *edgeMap[triangleIDs[i * 3 + 1]], *edgeMap[triangleIDs[i * 3 + 2]]);
      }
    }
  }
}

void CommunicateMesh::broadcastSendMesh(
    const mesh::Mesh &mesh)
{
  TRACE(mesh.getName());
  int dim = mesh.getDimensions();

  int numberOfVertices = mesh.vertices().size();
  _communication->broadcast(numberOfVertices);
  if (numberOfVertices > 0) {
    double coords[numberOfVertices * dim];
    int    globalIDs[numberOfVertices];
    for (int i = 0; i < numberOfVertices; i++) {
      for (int d = 0; d < dim; d++) {
        coords[i * dim + d] = mesh.vertices()[i].getCoords()[d];
      }
      globalIDs[i] = mesh.vertices()[i].getGlobalIndex();
    }
    _communication->broadcast(coords, numberOfVertices * dim);
    _communication->broadcast(globalIDs, numberOfVertices);
  }

  int numberOfEdges = mesh.edges().size();
  _communication->broadcast(numberOfEdges);
  if (numberOfEdges > 0) {
    //we need to send the vertexIDs first such that the right edges can be created later
    //contrary to the normal sendMesh, this variant must also work for adding delta meshes
    int vertexIDs[numberOfVertices];
    for (int i = 0; i < numberOfVertices; i++) {
      vertexIDs[i] = mesh.vertices()[i].getID();
    }
    _communication->broadcast(vertexIDs, numberOfVertices);

    int edgeIDs[numberOfEdges * 2];
    for (int i = 0; i < numberOfEdges; i++) {
      edgeIDs[i * 2]     = mesh.edges()[i].vertex(0).getID();
      edgeIDs[i * 2 + 1] = mesh.edges()[i].vertex(1).getID();
    }
    _communication->broadcast(edgeIDs, numberOfEdges * 2);
  }

  if (dim == 3) {
    int numberOfTriangles = mesh.triangles().size();
    _communication->broadcast(numberOfTriangles);
    if (numberOfTriangles > 0) {
      //we need to send the edgeIDs first such that the right edges can be created later
      //contrary to the normal sendMesh, this variant must also work for adding delta meshes
      int edgeIDs[numberOfEdges];
      for (int i = 0; i < numberOfEdges; i++) {
        edgeIDs[i] = mesh.edges()[i].getID();
      }
      _communication->broadcast(edgeIDs, numberOfEdges);

      int triangleIDs[numberOfTriangles * 3];
      for (int i = 0; i < numberOfTriangles; i++) {
        triangleIDs[i * 3]     = mesh.triangles()[i].edge(0).getID();
        triangleIDs[i * 3 + 1] = mesh.triangles()[i].edge(1).getID();
        triangleIDs[i * 3 + 2] = mesh.triangles()[i].edge(2).getID();
      }
      _communication->broadcast(triangleIDs, numberOfTriangles * 3);
    }
  }
}

void CommunicateMesh::broadcastReceiveMesh(
    mesh::Mesh &mesh)
{
  TRACE(mesh.getName());
  int dim             = mesh.getDimensions();
  int rankBroadcaster = 0;

  std::vector<mesh::Vertex *>   vertices;
  std::map<int, mesh::Vertex *> vertexMap;
  int                           numberOfVertices = 0;
  _communication->broadcast(numberOfVertices, rankBroadcaster);

  if (numberOfVertices > 0) {
    double vertexCoords[numberOfVertices * dim];
    int    globalIDs[numberOfVertices];
    _communication->broadcast(vertexCoords, numberOfVertices * dim, rankBroadcaster);
    _communication->broadcast(globalIDs, numberOfVertices, rankBroadcaster);
    for (int i = 0; i < numberOfVertices; i++) {
      Eigen::VectorXd coords(dim);
      for (int d = 0; d < dim; d++) {
        coords[d] = vertexCoords[i * dim + d];
      }
      mesh::Vertex &v = mesh.createVertex(coords);
      assertion(v.getID() >= 0, v.getID());
      v.setGlobalIndex(globalIDs[i]);
      vertices.push_back(&v);
    }
  }

  int                       numberOfEdges = 0;
  std::vector<mesh::Edge *> edges;
  _communication->broadcast(numberOfEdges, rankBroadcaster);
  if (numberOfEdges > 0) {
    int vertexIDs[numberOfVertices];
    _communication->broadcast(vertexIDs, numberOfVertices, rankBroadcaster);
    for (int i = 0; i < numberOfVertices; i++) {
      vertexMap[vertexIDs[i]] = vertices[i];
    }

    int edgeIDs[numberOfEdges * 2];
    _communication->broadcast(edgeIDs, numberOfEdges * 2, rankBroadcaster);
    for (int i = 0; i < numberOfEdges; i++) {
      assertion(vertexMap.find(edgeIDs[i * 2]) != vertexMap.end());
      assertion(vertexMap.find(edgeIDs[i * 2 + 1]) != vertexMap.end());
      assertion(edgeIDs[i * 2] != edgeIDs[i * 2 + 1]);
      mesh::Edge &e = mesh.createEdge(*vertexMap[edgeIDs[i * 2]], *vertexMap[edgeIDs[i * 2 + 1]]);
      edges.push_back(&e);
    }
  }

  if (dim == 3) {
    int numberOfTriangles = 0;
    _communication->broadcast(numberOfTriangles, rankBroadcaster);
    if (numberOfTriangles > 0) {
      assertion((edges.size() > 0) || (numberOfTriangles == 0));
      int edgeIDs[numberOfEdges];
      _communication->broadcast(edgeIDs, numberOfEdges, rankBroadcaster);
      std::map<int, mesh::Edge *> edgeMap;
      for (int i = 0; i < numberOfEdges; i++) {
        edgeMap[edgeIDs[i]] = edges[i];
      }

      int triangleIDs[numberOfTriangles * 3];
      _communication->broadcast(triangleIDs, numberOfTriangles * 3, rankBroadcaster);

      for (int i = 0; i < numberOfTriangles; i++) {
        assertion(edgeMap.find(triangleIDs[i * 3]) != edgeMap.end());
        assertion(edgeMap.find(triangleIDs[i * 3 + 1]) != edgeMap.end());
        assertion(edgeMap.find(triangleIDs[i * 3 + 2]) != edgeMap.end());
        assertion(triangleIDs[i * 3] != triangleIDs[i * 3 + 1]);
        assertion(triangleIDs[i * 3 + 1] != triangleIDs[i * 3 + 2]);
        assertion(triangleIDs[i * 3 + 2] != triangleIDs[i * 3]);
        mesh.createTriangle(*edgeMap[triangleIDs[i * 3]], *edgeMap[triangleIDs[i * 3 + 1]], *edgeMap[triangleIDs[i * 3 + 2]]);
      }
    }
  }
}

void CommunicateMesh::sendBoundingBox(
    const mesh::Mesh::BoundingBox &bb,
    int                            rankReceiver)
{
  TRACE(rankReceiver);
  int dim = bb.size();
  for (int d = 0; d < dim; d++) {
    _communication->send(bb[d].first, rankReceiver);
    _communication->send(bb[d].second, rankReceiver);
  }
}

void CommunicateMesh::receiveBoundingBox(
    mesh::Mesh::BoundingBox &bb,
    int                      rankSender)
{
  TRACE(rankSender);
  int dim = bb.size();
  for (int d = 0; d < dim; d++) {
    _communication->receive(bb[d].first, rankSender);
    _communication->receive(bb[d].second, rankSender);
  }
}
}
} // namespace precice, com
