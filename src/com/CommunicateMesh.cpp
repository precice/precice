#include <Eigen/Core>
#include <algorithm>
#include <boost/container/flat_map.hpp>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "CommunicateMesh.hpp"
#include "Communication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace com {
CommunicateMesh::CommunicateMesh(
    com::PtrCommunication communication)
    : _communication(std::move(communication))
{
}

void CommunicateMesh::sendMesh(
    const mesh::Mesh &mesh,
    int               rankReceiver)
{
  PRECICE_TRACE(mesh.getName(), rankReceiver);
  const int dim = mesh.getDimensions();

  const auto &meshVertices     = mesh.vertices();
  const int   numberOfVertices = meshVertices.size();
  _communication->send(numberOfVertices, rankReceiver);

  if (mesh.vertices().empty()) {
    return;
  }

  {
    std::vector<double> coords(static_cast<size_t>(numberOfVertices) * dim);
    std::vector<int>    globalIDs(numberOfVertices);
    for (size_t i = 0; i < static_cast<size_t>(numberOfVertices); ++i) {
      std::copy_n(meshVertices[i].rawCoords().begin(), dim, &coords[i * dim]);
      globalIDs[i] = meshVertices[i].getGlobalIndex();
    }
    _communication->sendRange(coords, rankReceiver);
    _communication->sendRange(globalIDs, rankReceiver);
  }

  _communication->send(mesh.hasConnectivity(), rankReceiver);
  if (!mesh.hasConnectivity()) {
    PRECICE_DEBUG("No connectivity to send");
    return;
  }

  // We need to send the vertexIDs first. This is required as the receiver will
  // end up with other vertexIDs after creating vertices.
  std::vector<int> vertexIDs(numberOfVertices);
  for (int i = 0; i < numberOfVertices; ++i) {
    vertexIDs[i] = meshVertices[i].getID();
  }
  _communication->sendRange(vertexIDs, rankReceiver);

  // Send Edges
  const int numberOfEdges = mesh.edges().size();
  PRECICE_DEBUG("Number of edges to send: {}", numberOfEdges);
  _communication->send(numberOfEdges, rankReceiver);
  if (mesh.hasEdges()) {
    std::vector<int> edgeIDs(numberOfEdges * 2);
    for (int i = 0; i < numberOfEdges; i++) {
      edgeIDs[i * 2]     = mesh.edges()[i].vertex(0).getID();
      edgeIDs[i * 2 + 1] = mesh.edges()[i].vertex(1).getID();
    }
    _communication->sendRange(edgeIDs, rankReceiver);
  }

  // Send Triangles
  int numberOfTriangles = mesh.triangles().size();
  PRECICE_DEBUG("Number of triangles to send: {}", numberOfTriangles);
  _communication->send(numberOfTriangles, rankReceiver);
  if (mesh.hasTriangles()) {
    std::vector<int> triangleIDs(numberOfTriangles * 3);
    for (int i = 0; i < numberOfTriangles; ++i) {
      triangleIDs[i * 3]     = mesh.triangles()[i].vertex(0).getID();
      triangleIDs[i * 3 + 1] = mesh.triangles()[i].vertex(1).getID();
      triangleIDs[i * 3 + 2] = mesh.triangles()[i].vertex(2).getID();
    }
    _communication->sendRange(triangleIDs, rankReceiver);
  }

  // Send Tetrahedra
  int numberOfTetra = mesh.tetrahedra().size();
  PRECICE_DEBUG("Number of tetrahedra to send: {}", numberOfTetra);
  _communication->send(numberOfTetra, rankReceiver);

  if (mesh.hasTetrahedra()) {

    std::vector<int> tetraIDs(numberOfTetra * 4);
    for (int i = 0; i < numberOfTetra; ++i) {
      tetraIDs[i * 4]     = mesh.tetrahedra()[i].vertex(0).getID();
      tetraIDs[i * 4 + 1] = mesh.tetrahedra()[i].vertex(1).getID();
      tetraIDs[i * 4 + 2] = mesh.tetrahedra()[i].vertex(2).getID();
      tetraIDs[i * 4 + 3] = mesh.tetrahedra()[i].vertex(3).getID();
    }
    _communication->sendRange(tetraIDs, rankReceiver);
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
  PRECICE_DEBUG("Number of vertices to receive: {}", numberOfVertices);

  if (numberOfVertices == 0) {
    return;
  }

  std::vector<mesh::Vertex *> vertices;
  vertices.reserve(numberOfVertices);
  {
    std::vector<double> vertexCoords = _communication->receiveRange(rankSender, AsVectorTag<double>{});
    std::vector<int>    globalIDs    = _communication->receiveRange(rankSender, AsVectorTag<int>{});
    Eigen::VectorXd     coords(dim);
    for (int i = 0; i < numberOfVertices; i++) {
      for (int d = 0; d < dim; d++) {
        coords[d] = vertexCoords[i * dim + d];
      }
      mesh::Vertex &v = mesh.createVertex(coords);
      PRECICE_ASSERT(v.getID() >= 0, v.getID());
      v.setGlobalIndex(globalIDs[i]);
      vertices.push_back(&v);
    }
  }

  bool hasConnectivity{false};
  _communication->receive(hasConnectivity, rankSender);
  if (!hasConnectivity) {
    PRECICE_DEBUG("No connectivity to receive");
    return;
  }

  // We need to receive the vertexIDs first. This is required as the vertices
  // created above have different vertexIDs as the original mesh. We need a mapping
  // from original to new vertexids.
  boost::container::flat_map<int, mesh::Vertex *> vertexMap;
  vertexMap.reserve(numberOfVertices);
  const std::vector<int> vertexIDs = _communication->receiveRange(rankSender, AsVectorTag<int>{});
  for (int i = 0; i < numberOfVertices; ++i) {
    vertexMap[vertexIDs[i]] = vertices[i];
  }

  // Receive Edges
  int numberOfEdges = 0;
  _communication->receive(numberOfEdges, rankSender);
  PRECICE_DEBUG("Number of edges to receive: {}", numberOfEdges);
  if (numberOfEdges > 0) {
    std::vector<int> edgeIDs = _communication->receiveRange(rankSender, AsVectorTag<int>{});
    for (int i = 0; i < numberOfEdges; i++) {
      PRECICE_ASSERT(vertexMap.count((edgeIDs[i * 2])) == 1);
      PRECICE_ASSERT(vertexMap.count(edgeIDs[i * 2 + 1]) == 1);
      PRECICE_ASSERT(edgeIDs[i * 2] != edgeIDs[i * 2 + 1]);
      mesh.createEdge(*vertexMap[edgeIDs[i * 2]], *vertexMap[edgeIDs[i * 2 + 1]]);
    }
  }

  // Receveive Triangles
  int numberOfTriangles = 0;
  _communication->receive(numberOfTriangles, rankSender);
  PRECICE_DEBUG("Number of triangles to receive: {}", numberOfTriangles);
  if (numberOfTriangles > 0) {
    std::vector<int> triangleIDs = _communication->receiveRange(rankSender, AsVectorTag<int>{});
    PRECICE_ASSERT(triangleIDs.size() == numberOfTriangles * 3);

    for (int i = 0; i < numberOfTriangles; i++) {
      PRECICE_ASSERT(vertexMap.count(triangleIDs[i * 3]) == 1);
      PRECICE_ASSERT(vertexMap.count(triangleIDs[i * 3 + 1]) == 1);
      PRECICE_ASSERT(vertexMap.count(triangleIDs[i * 3 + 2]) == 1);

      mesh.createTriangle(*vertexMap[triangleIDs[i * 3]], *vertexMap[triangleIDs[i * 3 + 1]], *vertexMap[triangleIDs[i * 3 + 2]]);
    }
  }

  // Receveive Tetrahedra

  int numberofTetra = 0;
  _communication->receive(numberofTetra, rankSender);
  PRECICE_DEBUG("Number of tetrahedra to receive: {}", numberofTetra);

  if (numberofTetra > 0) {
    std::vector<int> tetraIDs = _communication->receiveRange(rankSender, AsVectorTag<int>{});
    PRECICE_ASSERT(tetraIDs.size() == numberofTetra * 4);

    for (int i = 0; i < numberofTetra; i++) {
      PRECICE_ASSERT(vertexMap.count(tetraIDs[i * 4]) == 1);
      PRECICE_ASSERT(vertexMap.count(tetraIDs[i * 4 + 1]) == 1);
      PRECICE_ASSERT(vertexMap.count(tetraIDs[i * 4 + 2]) == 1);
      PRECICE_ASSERT(vertexMap.count(tetraIDs[i * 4 + 3]) == 1);

      mesh.createTetrahedron(*vertexMap[tetraIDs[i * 4]], *vertexMap[tetraIDs[i * 4 + 1]], *vertexMap[tetraIDs[i * 4 + 2]], *vertexMap[tetraIDs[i * 4 + 3]]);
    }
  }
}

void CommunicateMesh::broadcastSendMesh(const mesh::Mesh &mesh)
{
  PRECICE_TRACE(mesh.getName());
  int dim = mesh.getDimensions();

  const auto &meshVertices     = mesh.vertices();
  const int   numberOfVertices = meshVertices.size();
  _communication->broadcast(numberOfVertices);
  if (numberOfVertices > 0) {
    std::vector<double> coords(static_cast<size_t>(numberOfVertices) * dim);
    std::vector<int>    globalIDs(numberOfVertices);
    for (size_t i = 0; i < static_cast<size_t>(numberOfVertices); ++i) {
      std::copy_n(meshVertices[i].rawCoords().begin(), dim, &coords[i * dim]);
      globalIDs[i] = meshVertices[i].getGlobalIndex();
    }
    _communication->broadcast(coords);
    _communication->broadcast(globalIDs);
  }

  _communication->broadcast(mesh.hasConnectivity());
  if (!mesh.hasConnectivity()) {
    return;
  }

  // We need to send the vertexIDs first. This is required as the receiver will
  // end up with other vertexIDs after creating vertices.
  std::vector<int> vertexIDs(numberOfVertices);
  for (int i = 0; i < numberOfVertices; i++) {
    vertexIDs[i] = meshVertices[i].getID();
  }
  _communication->broadcast(vertexIDs);

  // Send Edges
  int numberOfEdges = mesh.edges().size();
  _communication->broadcast(numberOfEdges);
  if (numberOfEdges > 0) {
    std::vector<int> edgeIDs(numberOfEdges * 2);
    const auto &     meshEdges = mesh.edges();
    for (int i = 0; i < numberOfEdges; i++) {
      edgeIDs[i * 2]     = meshEdges[i].vertex(0).getID();
      edgeIDs[i * 2 + 1] = meshEdges[i].vertex(1).getID();
    }
    _communication->broadcast(edgeIDs);
  }

  // Send Triangles
  int numberOfTriangles = mesh.triangles().size();
  _communication->broadcast(numberOfTriangles);
  if (numberOfTriangles > 0) {
    std::vector<int> triangleIDs(numberOfTriangles * 3);
    const auto &     meshTriangles = mesh.triangles();
    for (int i = 0; i < numberOfTriangles; i++) {
      triangleIDs[i * 3]     = meshTriangles[i].vertex(0).getID();
      triangleIDs[i * 3 + 1] = meshTriangles[i].vertex(1).getID();
      triangleIDs[i * 3 + 2] = meshTriangles[i].vertex(2).getID();
    }
    _communication->broadcast(triangleIDs);
  }

  // Send Tetrahedra
  int numberOfTetra = mesh.tetrahedra().size();
  _communication->broadcast(numberOfTetra);

  if (numberOfTetra > 0) {
    std::vector<int> tetraIDs(numberOfTetra * 4);
    const auto &     meshTetrahedra = mesh.tetrahedra();
    for (int i = 0; i < numberOfTetra; i++) {
      tetraIDs[i * 4]     = meshTetrahedra[i].vertex(0).getID();
      tetraIDs[i * 4 + 1] = meshTetrahedra[i].vertex(1).getID();
      tetraIDs[i * 4 + 2] = meshTetrahedra[i].vertex(2).getID();
      tetraIDs[i * 4 + 3] = meshTetrahedra[i].vertex(2).getID();
    }
    _communication->broadcast(tetraIDs);
  }
}

void CommunicateMesh::broadcastReceiveMesh(
    mesh::Mesh &mesh)
{
  PRECICE_TRACE(mesh.getName());
  int  dim             = mesh.getDimensions();
  Rank rankBroadcaster = 0;

  std::vector<mesh::Vertex *> vertices;
  int                         numberOfVertices = 0;
  _communication->broadcast(numberOfVertices, rankBroadcaster);

  if (numberOfVertices > 0) {
    std::vector<double> vertexCoords;
    std::vector<int>    globalIDs;
    _communication->broadcast(vertexCoords, rankBroadcaster);
    _communication->broadcast(globalIDs, rankBroadcaster);
    Eigen::VectorXd coords(dim);
    for (int i = 0; i < numberOfVertices; i++) {
      for (int d = 0; d < dim; d++) {
        coords[d] = vertexCoords[i * dim + d];
      }
      mesh::Vertex &v = mesh.createVertex(coords);
      PRECICE_ASSERT(v.getID() >= 0, v.getID());
      v.setGlobalIndex(globalIDs[i]);
      vertices.push_back(&v);
    }
  }

  bool hasConnectivity{false};
  _communication->broadcast(hasConnectivity, rankBroadcaster);
  if (!hasConnectivity) {
    return;
  }

  // We need to receive the vertexIDs first. This is required as the vertices
  // created above have different vertexIDs as the original mesh. We need a mapping
  // from original to new vertexids.
  std::vector<int> vertexIDs;
  _communication->broadcast(vertexIDs, rankBroadcaster);
  boost::container::flat_map<VertexID, mesh::Vertex *> vertexMap;
  vertexMap.reserve(vertexIDs.size());
  for (int i = 0; i < numberOfVertices; i++) {
    vertexMap[vertexIDs[i]] = vertices[i];
  }

  // Receive Edges
  int numberOfEdges = 0;
  _communication->broadcast(numberOfEdges, rankBroadcaster);
  if (numberOfEdges > 0) {
    std::vector<int> edgeIDs;
    _communication->broadcast(edgeIDs, rankBroadcaster);
    for (int i = 0; i < numberOfEdges; i++) {
      PRECICE_ASSERT(vertexMap.count(edgeIDs[i * 2]) == 1);
      PRECICE_ASSERT(vertexMap.count(edgeIDs[i * 2 + 1]) == 1);
      PRECICE_ASSERT(edgeIDs[i * 2] != edgeIDs[i * 2 + 1]);
      mesh.createEdge(*vertexMap[edgeIDs[i * 2]], *vertexMap[edgeIDs[i * 2 + 1]]);
    }
  }

  // Receive Triangles
  int numberOfTriangles = 0;
  _communication->broadcast(numberOfTriangles, rankBroadcaster);
  if (numberOfTriangles > 0) {
    std::vector<int> triangleIDs;
    _communication->broadcast(triangleIDs, rankBroadcaster);
    for (int i = 0; i < numberOfTriangles; i++) {
      PRECICE_ASSERT(vertexMap.count(triangleIDs[i * 3]) == 1);
      PRECICE_ASSERT(vertexMap.count(triangleIDs[i * 3 + 1]) == 1);
      PRECICE_ASSERT(vertexMap.count(triangleIDs[i * 3 + 2]) == 1);
      PRECICE_ASSERT(triangleIDs[i * 3] != triangleIDs[i * 3 + 1]);
      PRECICE_ASSERT(triangleIDs[i * 3 + 1] != triangleIDs[i * 3 + 2]);
      PRECICE_ASSERT(triangleIDs[i * 3 + 2] != triangleIDs[i * 3]);
      mesh.createTriangle(*vertexMap[triangleIDs[i * 3]], *vertexMap[triangleIDs[i * 3 + 1]], *vertexMap[triangleIDs[i * 3 + 2]]);
    }
  }

  // Receive Tetrahedra
  int numberOfTetra = 0;
  _communication->broadcast(numberOfTetra, rankBroadcaster);

  if (numberOfTetra > 0) {
    std::vector<int> tetraIDs;
    _communication->broadcast(tetraIDs, rankBroadcaster);
    for (int i = 0; i < numberOfTetra; i++) {
      PRECICE_ASSERT(vertexMap.count(tetraIDs[i * 4]) == 1);
      PRECICE_ASSERT(vertexMap.count(tetraIDs[i * 4 + 1]) == 1);
      PRECICE_ASSERT(vertexMap.count(tetraIDs[i * 4 + 2]) == 1);
      PRECICE_ASSERT(vertexMap.count(tetraIDs[i * 4 + 3]) == 1);

      PRECICE_ASSERT(tetraIDs[i * 4] != tetraIDs[i * 4 + 1]);
      PRECICE_ASSERT(tetraIDs[i * 4 + 1] != tetraIDs[i * 4 + 2]);
      PRECICE_ASSERT(tetraIDs[i * 4 + 2] != tetraIDs[i * 4]);
      PRECICE_ASSERT(tetraIDs[i * 4 + 3] != tetraIDs[i * 4]);
      PRECICE_ASSERT(tetraIDs[i * 4 + 3] != tetraIDs[i * 4 + 1]);
      PRECICE_ASSERT(tetraIDs[i * 4 + 3] != tetraIDs[i * 4 + 2]);
      mesh.createTetrahedron(*vertexMap[tetraIDs[i * 4]], *vertexMap[tetraIDs[i * 4 + 1]], *vertexMap[tetraIDs[i * 4 + 2]], *vertexMap[tetraIDs[i * 4 + 3]]);
    }
  }
}

} // namespace com
} // namespace precice
