#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "com/CommunicateMesh.hpp"
#include "com/SharedPointer.hpp"
#include "m2n/M2N.hpp"
#include "mesh/Mesh.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace mesh {
class Edge;
class Triangle;
class Vertex;
} // namespace mesh
} // namespace precice

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(VertexEdgeMesh)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  auto m2n = context.connectMasters("A", "B");

  for (int dim = 2; dim <= 3; dim++) {
    mesh::Mesh    sendMesh("Sent Mesh", dim, false, testing::nextMeshID());
    mesh::Vertex &v0 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 0));
    mesh::Vertex &v1 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
    mesh::Vertex &v2 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 2));
    mesh::Edge &  e0 = sendMesh.createEdge(v0, v1);
    mesh::Edge &  e1 = sendMesh.createEdge(v1, v2);
    mesh::Edge &  e2 = sendMesh.createEdge(v2, v0);

    CommunicateMesh comMesh(m2n->getMasterCommunication());

    if (context.isNamed("A")) {
      comMesh.sendMesh(sendMesh, 0);
    } else {
      // receiveMesh can also deal with delta meshes
      mesh::Mesh recvMesh("Received Mesh", dim, false, testing::nextMeshID());
      recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
      comMesh.receiveMesh(recvMesh, 0);
      BOOST_TEST(recvMesh.vertices().size() == 4);
      BOOST_TEST(testing::equals(recvMesh.vertices().at(0).getCoords(), Eigen::VectorXd::Constant(dim, 9)));
      BOOST_TEST(recvMesh.vertices().at(1) == v0);
      BOOST_TEST(recvMesh.vertices().at(2) == v1);
      BOOST_TEST(recvMesh.vertices().at(3) == v2);
      BOOST_TEST(recvMesh.edges().at(0) == e0);
      BOOST_TEST(recvMesh.edges().at(1) == e1);
      BOOST_TEST(recvMesh.edges().at(2) == e2);
    }
  }
}

BOOST_AUTO_TEST_CASE(VertexEdgeTriangleMesh)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), Require::Events);
  auto m2n = context.connectMasters("A", "B");

  int             dim = 3;
  mesh::Mesh      sendMesh("Sent Mesh", dim, false, testing::nextMeshID());
  mesh::Vertex &  v0 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 0));
  mesh::Vertex &  v1 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
  mesh::Vertex &  v2 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 2));
  mesh::Edge &    e0 = sendMesh.createEdge(v0, v1);
  mesh::Edge &    e1 = sendMesh.createEdge(v1, v2);
  mesh::Edge &    e2 = sendMesh.createEdge(v2, v0);
  mesh::Triangle &t0 = sendMesh.createTriangle(e0, e1, e2);

  // Create mesh communicator
  CommunicateMesh comMesh(m2n->getMasterCommunication());

  if (context.isNamed("A")) {
    comMesh.sendMesh(sendMesh, 0);
  } else {
    mesh::Mesh recvMesh("Received Mesh", dim, false, testing::nextMeshID());
    // receiveMesh can also deal with delta meshes
    recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
    comMesh.receiveMesh(recvMesh, 0);
    BOOST_TEST(recvMesh.vertices().size() == 4);
    BOOST_TEST(testing::equals(recvMesh.vertices().at(0).getCoords(), Eigen::VectorXd::Constant(dim, 9)));
    BOOST_TEST(recvMesh.vertices().at(1) == v0);
    BOOST_TEST(recvMesh.vertices().at(2) == v1);
    BOOST_TEST(recvMesh.vertices().at(3) == v2);
    BOOST_TEST(recvMesh.edges().at(0) == e0);
    BOOST_TEST(recvMesh.edges().at(1) == e1);
    BOOST_TEST(recvMesh.edges().at(2) == e2);

    BOOST_TEST(recvMesh.triangles().at(0) == t0);
  }
}

BOOST_AUTO_TEST_CASE(BroadcastVertexEdgeTriangleMesh)
{
  PRECICE_TEST(""_on(2_ranks).setupMasterSlaves(), Require::Events);

  int             dim = 3;
  mesh::Mesh      sendMesh("Sent Mesh", dim, false, testing::nextMeshID());
  mesh::Vertex &  v0 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 0));
  mesh::Vertex &  v1 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
  mesh::Vertex &  v2 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 2));
  mesh::Edge &    e0 = sendMesh.createEdge(v0, v1);
  mesh::Edge &    e1 = sendMesh.createEdge(v1, v2);
  mesh::Edge &    e2 = sendMesh.createEdge(v2, v0);
  mesh::Triangle &t0 = sendMesh.createTriangle(e0, e1, e2);

  // Create mesh communicator
  CommunicateMesh comMesh(precice::utils::MasterSlave::_communication);

  if (context.isMaster()) {
    comMesh.broadcastSendMesh(sendMesh);
  } else {
    mesh::Mesh recvMesh("Received Mesh", dim, false, testing::nextMeshID());
    // receiveMesh can also deal with delta meshes
    recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
    comMesh.broadcastReceiveMesh(recvMesh);
    BOOST_TEST(recvMesh.vertices().size() == 4);
    BOOST_TEST(testing::equals(recvMesh.vertices().at(0).getCoords(), Eigen::VectorXd::Constant(dim, 9)));
    BOOST_TEST(recvMesh.vertices().at(1) == v0);
    BOOST_TEST(recvMesh.vertices().at(2) == v1);
    BOOST_TEST(recvMesh.vertices().at(3) == v2);
    BOOST_TEST(recvMesh.edges().at(0) == e0);
    BOOST_TEST(recvMesh.edges().at(1) == e1);
    BOOST_TEST(recvMesh.edges().at(2) == e2);
    BOOST_TEST(recvMesh.triangles().at(0) == t0);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Mesh
BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI
