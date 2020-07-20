#include "CommunicateMesh.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <boost/container/flat_map.hpp>
#include <functional>
#include <map>
#include <memory>
#include <ostream>
#include <stddef.h>
#include <vector>
#include "Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace com {
CommunicateMesh::CommunicateMesh(
    com::PtrCommunication communication)
    : _communication(communication)
{
}

void CommunicateMesh::sendMesh(
    const mesh::Mesh &mesh,
    int               rankReceiver)
{
  PRECICE_TRACE(mesh.getName(), rankReceiver);
  int dim = mesh.getDimensions();

  int numberOfVertices = mesh.vertices().size();
  _communication->send(numberOfVertices, rankReceiver);
  if (not mesh.vertices().empty()) {
    std::vector<double> coords(static_cast<size_t>(numberOfVertices) * dim);
    std::vector<int>    globalIDs(numberOfVertices);
    for (int i = 0; i < numberOfVertices; i++) {
      for (int d = 0; d < dim; d++) {
        coords[i * dim + d] = mesh.vertices()[i].getCoords()[d];
      }
      globalIDs[i] = mesh.vertices()[i].getGlobalIndex();
    }
    _communication->send(coords, rankReceiver);
    _communication->send(globalIDs, rankReceiver);
  }

  int numberOfEdges = mesh.edges().size();
  _communication->send(numberOfEdges, rankReceiver);
  if (not mesh.edges().empty()) {
    //we need to send the vertexIDs first such that the right edges can be created later
    //contrary to the normal sendMesh, this variant must also work for adding delta meshes
    std::vector<int> vertexIDs(numberOfVertices);
    for (int i = 0; i < numberOfVertices; i++) {
      vertexIDs[i] = mesh.vertices()[i].getID();
    }
    _communication->send(vertexIDs, rankReceiver);

    std::vector<int> edgeIDs(numberOfEdges * 2);
    for (int i = 0; i < numberOfEdges; i++) {
      edgeIDs[i * 2]     = mesh.edges()[i].vertex(0).getID();
      edgeIDs[i * 2 + 1] = mesh.edges()[i].vertex(1).getID();
    }
    _communication->send(edgeIDs, rankReceiver);
  }

  if (dim == 3) {
    int numberOfTriangles = mesh.triangles().size();
    _communication->send(numberOfTriangles, rankReceiver);
    if (not mesh.triangles().empty()) {
      //we need to send the edgeIDs first such that the right edges can be created later
      //contrary to the normal sendMesh, this variant must also work for adding delta meshes
      std::vector<int> edgeIDs(numberOfEdges);
      for (int i = 0; i < numberOfEdges; i++) {
        edgeIDs[i] = mesh.edges()[i].getID();
      }
      _communication->send(edgeIDs, rankReceiver);

      std::vector<int> triangleIDs(numberOfTriangles * 3);
      for (int i = 0; i < numberOfTriangles; i++) {
        triangleIDs[i * 3]     = mesh.triangles()[i].edge(0).getID();
        triangleIDs[i * 3 + 1] = mesh.triangles()[i].edge(1).getID();
        triangleIDs[i * 3 + 2] = mesh.triangles()[i].edge(2).getID();
      }
      _communication->send(triangleIDs, rankReceiver);
    }
  }
}

void CommunicateMesh::receiveMesh(
    mesh::Mesh &mesh,
    int         rankSender)
{
  PRECICE_TRACE(mesh.getName(), rankSender);
  int dim = mesh.getDimensions();

  int numberOfVertices = 0;
  _communication->receive(numberOfVertices, rankSender);
  PRECICE_DEBUG("Number of vertices to receive: " << numberOfVertices);

  std::vector<mesh::Vertex *> vertices;
  vertices.reserve(numberOfVertices);
  if (numberOfVertices > 0) {
    std::vector<double> vertexCoords;
    std::vector<int>    globalIDs;
    _communication->receive(vertexCoords, rankSender);
    _communication->receive(globalIDs, rankSender);
    for (int i = 0; i < numberOfVertices; i++) {
      Eigen::VectorXd coords(dim);
      for (int d = 0; d < dim; d++) {
        coords[d] = vertexCoords[i * dim + d];
      }
      mesh::Vertex &v = mesh.createVertex(coords);
      PRECICE_ASSERT(v.getID() >= 0, v.getID());
      v.setGlobalIndex(globalIDs[i]);
      vertices.push_back(&v);
    }
  }

  int numberOfEdges = 0;
  _communication->receive(numberOfEdges, rankSender);
  PRECICE_DEBUG("Number of edges to receive: " << numberOfEdges);

  boost::container::flat_map<int, mesh::Vertex *> vertexMap;
  vertexMap.reserve(numberOfVertices);
  std::vector<mesh::Edge *> edges;
  if (numberOfEdges > 0) {
    std::vector<int> vertexIDs;
    _communication->receive(vertexIDs, rankSender);
    for (int i = 0; i < numberOfVertices; i++) {
      vertexMap[vertexIDs[i]] = vertices[i];
    }

    std::vector<int> edgeIDs;
    _communication->receive(edgeIDs, rankSender);
    for (int i = 0; i < numberOfEdges; i++) {
      PRECICE_ASSERT(vertexMap.count((edgeIDs[i * 2])) == 1);
      PRECICE_ASSERT(vertexMap.count(edgeIDs[i * 2 + 1]) == 1);
      PRECICE_ASSERT(edgeIDs[i * 2] != edgeIDs[i * 2 + 1]);
      mesh::Edge &e = mesh.createEdge(*vertexMap[edgeIDs[i * 2]], *vertexMap[edgeIDs[i * 2 + 1]]);
      edges.push_back(&e);
    }
  }

  if (dim == 3) {
    int numberOfTriangles = 0;
    _communication->receive(numberOfTriangles, rankSender);
    PRECICE_DEBUG("Number of Triangles to receive: " << numberOfTriangles);
    PRECICE_DEBUG("Number of Edges: " << edges.size());
    if (numberOfTriangles > 0) {
      PRECICE_ASSERT((edges.size() > 0) || (numberOfTriangles == 0));
      std::vector<int> edgeIDs;
      _communication->receive(edgeIDs, rankSender);
      boost::container::flat_map<int, mesh::Edge *> edgeMap;
      edgeMap.reserve(numberOfEdges);
      for (int i = 0; i < numberOfEdges; i++) {
        edgeMap[edgeIDs[i]] = edges[i];
      }

      std::vector<int> triangleIDs;
      _communication->receive(triangleIDs, rankSender);

      for (int i = 0; i < numberOfTriangles; i++) {
        PRECICE_ASSERT(edgeMap.count(triangleIDs[i * 3]) == 1);
        PRECICE_ASSERT(edgeMap.count(triangleIDs[i * 3 + 1]) == 1);
        PRECICE_ASSERT(edgeMap.count(triangleIDs[i * 3 + 2]) == 1);
        PRECICE_ASSERT(triangleIDs[i * 3] != triangleIDs[i * 3 + 1]);
        PRECICE_ASSERT(triangleIDs[i * 3 + 1] != triangleIDs[i * 3 + 2]);
        PRECICE_ASSERT(triangleIDs[i * 3 + 2] != triangleIDs[i * 3]);
        mesh.createTriangle(*edgeMap[triangleIDs[i * 3]], *edgeMap[triangleIDs[i * 3 + 1]], *edgeMap[triangleIDs[i * 3 + 2]]);
      }
    }
  }
}

void CommunicateMesh::broadcastSendMesh(const mesh::Mesh &mesh)
{
  PRECICE_TRACE(mesh.getName());
  int dim = mesh.getDimensions();

  int numberOfVertices = mesh.vertices().size();
  _communication->broadcast(numberOfVertices);
  if (numberOfVertices > 0) {
    std::vector<double> coords(static_cast<size_t>(numberOfVertices) * dim);
    std::vector<int>    globalIDs(numberOfVertices);
    for (int i = 0; i < numberOfVertices; i++) {
      for (int d = 0; d < dim; d++) {
        coords[i * dim + d] = mesh.vertices()[i].getCoords()[d];
      }
      globalIDs[i] = mesh.vertices()[i].getGlobalIndex();
    }
    _communication->broadcast(coords);
    _communication->broadcast(globalIDs);
  }

  int numberOfEdges = mesh.edges().size();
  _communication->broadcast(numberOfEdges);
  if (numberOfEdges > 0) {
    //we need to send the vertexIDs first such that the right edges can be created later
    //contrary to the normal sendMesh, this variant must also work for adding delta meshes
    std::vector<int> vertexIDs(numberOfVertices);
    for (int i = 0; i < numberOfVertices; i++) {
      vertexIDs[i] = mesh.vertices()[i].getID();
    }
    _communication->broadcast(vertexIDs);

    std::vector<int> edgeIDs(numberOfEdges * 2);
    for (int i = 0; i < numberOfEdges; i++) {
      edgeIDs[i * 2]     = mesh.edges()[i].vertex(0).getID();
      edgeIDs[i * 2 + 1] = mesh.edges()[i].vertex(1).getID();
    }
    _communication->broadcast(edgeIDs);
  }

  if (dim == 3) {
    int numberOfTriangles = mesh.triangles().size();
    _communication->broadcast(numberOfTriangles);
    if (numberOfTriangles > 0) {
      //we need to send the edgeIDs first such that the right edges can be created later
      //contrary to the normal sendMesh, this variant must also work for adding delta meshes
      std::vector<int> edgeIDs(numberOfEdges);
      for (int i = 0; i < numberOfEdges; i++) {
        edgeIDs[i] = mesh.edges()[i].getID();
      }
      _communication->broadcast(edgeIDs);

      std::vector<int> triangleIDs(numberOfTriangles * 3);
      for (int i = 0; i < numberOfTriangles; i++) {
        triangleIDs[i * 3]     = mesh.triangles()[i].edge(0).getID();
        triangleIDs[i * 3 + 1] = mesh.triangles()[i].edge(1).getID();
        triangleIDs[i * 3 + 2] = mesh.triangles()[i].edge(2).getID();
      }
      _communication->broadcast(triangleIDs);
    }
  }
}

void CommunicateMesh::broadcastReceiveMesh(
    mesh::Mesh &mesh)
{
  PRECICE_TRACE(mesh.getName());
  int dim             = mesh.getDimensions();
  int rankBroadcaster = 0;

  std::vector<mesh::Vertex *>   vertices;
  std::map<int, mesh::Vertex *> vertexMap;
  int                           numberOfVertices = 0;
  _communication->broadcast(numberOfVertices, rankBroadcaster);

  if (numberOfVertices > 0) {
    std::vector<double> vertexCoords;
    std::vector<int>    globalIDs;
    _communication->broadcast(vertexCoords, rankBroadcaster);
    _communication->broadcast(globalIDs, rankBroadcaster);
    for (int i = 0; i < numberOfVertices; i++) {
      Eigen::VectorXd coords(dim);
      for (int d = 0; d < dim; d++) {
        coords[d] = vertexCoords[i * dim + d];
      }
      mesh::Vertex &v = mesh.createVertex(coords);
      PRECICE_ASSERT(v.getID() >= 0, v.getID());
      v.setGlobalIndex(globalIDs[i]);
      vertices.push_back(&v);
    }
  }

  int                       numberOfEdges = 0;
  std::vector<mesh::Edge *> edges;
  _communication->broadcast(numberOfEdges, rankBroadcaster);
  if (numberOfEdges > 0) {
    std::vector<int> vertexIDs;
    _communication->broadcast(vertexIDs, rankBroadcaster);
    for (int i = 0; i < numberOfVertices; i++) {
      vertexMap[vertexIDs[i]] = vertices[i];
    }

    std::vector<int> edgeIDs;
    _communication->broadcast(edgeIDs, rankBroadcaster);
    for (int i = 0; i < numberOfEdges; i++) {
      PRECICE_ASSERT(vertexMap.find(edgeIDs[i * 2]) != vertexMap.end());
      PRECICE_ASSERT(vertexMap.find(edgeIDs[i * 2 + 1]) != vertexMap.end());
      PRECICE_ASSERT(edgeIDs[i * 2] != edgeIDs[i * 2 + 1]);
      mesh::Edge &e = mesh.createEdge(*vertexMap[edgeIDs[i * 2]], *vertexMap[edgeIDs[i * 2 + 1]]);
      edges.push_back(&e);
    }
  }

  if (dim == 3) {
    int numberOfTriangles = 0;
    _communication->broadcast(numberOfTriangles, rankBroadcaster);
    if (numberOfTriangles > 0) {
      PRECICE_ASSERT((edges.size() > 0) || (numberOfTriangles == 0));
      std::vector<int> edgeIDs;
      _communication->broadcast(edgeIDs, rankBroadcaster);
      std::map<int, mesh::Edge *> edgeMap;
      for (int i = 0; i < numberOfEdges; i++) {
        edgeMap[edgeIDs[i]] = edges[i];
      }

      std::vector<int> triangleIDs;
      _communication->broadcast(triangleIDs, rankBroadcaster);

      for (int i = 0; i < numberOfTriangles; i++) {
        PRECICE_ASSERT(edgeMap.find(triangleIDs[i * 3]) != edgeMap.end());
        PRECICE_ASSERT(edgeMap.find(triangleIDs[i * 3 + 1]) != edgeMap.end());
        PRECICE_ASSERT(edgeMap.find(triangleIDs[i * 3 + 2]) != edgeMap.end());
        PRECICE_ASSERT(triangleIDs[i * 3] != triangleIDs[i * 3 + 1]);
        PRECICE_ASSERT(triangleIDs[i * 3 + 1] != triangleIDs[i * 3 + 2]);
        PRECICE_ASSERT(triangleIDs[i * 3 + 2] != triangleIDs[i * 3]);
        mesh.createTriangle(*edgeMap[triangleIDs[i * 3]], *edgeMap[triangleIDs[i * 3 + 1]], *edgeMap[triangleIDs[i * 3 + 2]]);
      }
    }
  }
}

} // namespace com
} // namespace precice
