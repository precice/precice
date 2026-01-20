#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "com/Extra.hpp"
#include "com/SharedPointer.hpp"
#include "m2n/M2N.hpp"
#include "mesh/Mesh.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/IntraComm.hpp"

namespace precice::mesh {
class Edge;
class Triangle;
class Vertex;
} // namespace precice::mesh

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MeshTests)

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(VertexEdgeMesh)
{
  PRECICE_TEST();
  auto m2n = context.connectPrimaryRanks("A", "B");

  for (int dim = 2; dim <= 3; dim++) {
    mesh::Mesh    sendMesh("Sent Mesh", dim, testing::nextMeshID());
    mesh::Vertex &v0 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 0));
    mesh::Vertex &v1 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
    mesh::Vertex &v2 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 2));
    mesh::Edge   &e0 = sendMesh.createEdge(v0, v1);
    mesh::Edge   &e1 = sendMesh.createEdge(v1, v2);
    mesh::Edge   &e2 = sendMesh.createEdge(v2, v0);

    auto &comm = *m2n->getPrimaryRankCommunication();

    if (context.isNamed("A")) {
      com::sendMesh(comm, 0, sendMesh);
    } else {
      // receiveMesh can also deal with delta meshes
      mesh::Mesh recvMesh("Received Mesh", dim, testing::nextMeshID());
      recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
      com::receiveMesh(comm, 0, recvMesh);
      BOOST_TEST(recvMesh.nVertices() == 4);
      BOOST_TEST(testing::equals(recvMesh.vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 9)));
      BOOST_TEST(recvMesh.vertex(1) == v0);
      BOOST_TEST(recvMesh.vertex(2) == v1);
      BOOST_TEST(recvMesh.vertex(3) == v2);
      BOOST_TEST(recvMesh.edges().at(0) == e0);
      BOOST_TEST(recvMesh.edges().at(1) == e1);
      BOOST_TEST(recvMesh.edges().at(2) == e2);
    }
  }
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(VertexEdgeTriangleMesh)
{
  PRECICE_TEST();
  auto m2n = context.connectPrimaryRanks("A", "B");

  int             dim = 3;
  mesh::Mesh      sendMesh("Sent Mesh", dim, testing::nextMeshID());
  mesh::Vertex   &v0 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 0));
  mesh::Vertex   &v1 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
  mesh::Vertex   &v2 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 2));
  mesh::Edge     &e0 = sendMesh.createEdge(v0, v1);
  mesh::Edge     &e1 = sendMesh.createEdge(v1, v2);
  mesh::Edge     &e2 = sendMesh.createEdge(v2, v0);
  mesh::Triangle &t0 = sendMesh.createTriangle(e0, e1, e2);

  // Create mesh communicator
  auto &comm = *m2n->getPrimaryRankCommunication();

  if (context.isNamed("A")) {
    com::sendMesh(comm, 0, sendMesh);
  } else {
    mesh::Mesh recvMesh("Received Mesh", dim, testing::nextMeshID());
    // receiveMesh can also deal with delta meshes
    recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
    com::receiveMesh(comm, 0, recvMesh);
    BOOST_TEST(recvMesh.nVertices() == 4);
    BOOST_TEST(testing::equals(recvMesh.vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 9)));
    BOOST_TEST(recvMesh.vertex(1) == v0);
    BOOST_TEST(recvMesh.vertex(2) == v1);
    BOOST_TEST(recvMesh.vertex(3) == v2);
    BOOST_TEST(recvMesh.edges().at(0) == e0);
    BOOST_TEST(recvMesh.edges().at(1) == e1);
    BOOST_TEST(recvMesh.edges().at(2) == e2);

    BOOST_TEST(recvMesh.triangles().at(0) == t0);
  }
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastVertexEdgeTriangleMesh)
{
  PRECICE_TEST();

  int             dim = 3;
  mesh::Mesh      sendMesh("Sent Mesh", dim, testing::nextMeshID());
  mesh::Vertex   &v0 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 0));
  mesh::Vertex   &v1 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
  mesh::Vertex   &v2 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 2));
  mesh::Edge     &e0 = sendMesh.createEdge(v0, v1);
  mesh::Edge     &e1 = sendMesh.createEdge(v1, v2);
  mesh::Edge     &e2 = sendMesh.createEdge(v2, v0);
  mesh::Triangle &t0 = sendMesh.createTriangle(e0, e1, e2);

  // Create mesh communicator
  auto &comm = *precice::utils::IntraComm::getCommunication();

  if (context.isPrimary()) {
    com::broadcastSendMesh(comm, sendMesh);
  } else {
    mesh::Mesh recvMesh("Received Mesh", dim, testing::nextMeshID());
    // receiveMesh can also deal with delta meshes
    recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
    com::broadcastReceiveMesh(comm, recvMesh);
    BOOST_TEST(recvMesh.nVertices() == 4);
    BOOST_TEST(testing::equals(recvMesh.vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 9)));
    BOOST_TEST(recvMesh.vertex(1) == v0);
    BOOST_TEST(recvMesh.vertex(2) == v1);
    BOOST_TEST(recvMesh.vertex(3) == v2);
    BOOST_TEST(recvMesh.edges().at(0) == e0);
    BOOST_TEST(recvMesh.edges().at(1) == e1);
    BOOST_TEST(recvMesh.edges().at(2) == e2);
    BOOST_TEST(recvMesh.triangles().at(0) == t0);
  }
}

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), Require::Events)
BOOST_AUTO_TEST_CASE(OneTetraCommunication)
{
  PRECICE_TEST();
  auto m2n = context.connectPrimaryRanks("A", "B");

  int           dim = 3;
  mesh::Mesh    sendMesh("Sent Mesh", dim, testing::nextMeshID());
  mesh::Vertex &v0 = sendMesh.createVertex(Eigen::Vector3d{0.0, 0.0, 0.0});
  mesh::Vertex &v1 = sendMesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
  mesh::Vertex &v2 = sendMesh.createVertex(Eigen::Vector3d{0.0, 1.0, 0.0});
  mesh::Vertex &v3 = sendMesh.createVertex(Eigen::Vector3d{0.0, 0.0, 1.0});

  mesh::Tetrahedron &t0 = sendMesh.createTetrahedron(v0, v1, v2, v3);

  // Create mesh communicator
  auto &comm = *m2n->getPrimaryRankCommunication();

  if (context.isNamed("A")) {
    com::sendMesh(comm, 0, sendMesh);
  } else {
    mesh::Mesh recvMesh("Received Mesh", dim, testing::nextMeshID());
    // receiveMesh can also deal with delta meshes
    recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
    com::receiveMesh(comm, 0, recvMesh);
    BOOST_TEST(recvMesh.nVertices() == 5); // 4 + 1
    BOOST_TEST(testing::equals(recvMesh.vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 9)));
    BOOST_TEST(recvMesh.tetrahedra().size() == 1);
    BOOST_TEST(testing::equals(recvMesh.tetrahedra()[0].vertex(0).getCoords(), Eigen::Vector3d{0.0, 0.0, 0.0}));
    BOOST_TEST(recvMesh.tetrahedra()[0] == t0);
  }
}

PRECICE_TEST_SETUP(""_on(2_ranks).setupIntraComm(), Require::Events)
BOOST_AUTO_TEST_CASE(BroadcastTetra)
{
  PRECICE_TEST();

  int           dim = 3;
  mesh::Mesh    sendMesh("Sent Mesh", dim, testing::nextMeshID());
  mesh::Vertex &v0 = sendMesh.createVertex(Eigen::Vector3d{0.0, 0.0, 0.0});
  mesh::Vertex &v1 = sendMesh.createVertex(Eigen::Vector3d{1.0, 0.0, 0.0});
  mesh::Vertex &v2 = sendMesh.createVertex(Eigen::Vector3d{0.0, 1.0, 0.0});
  mesh::Vertex &v3 = sendMesh.createVertex(Eigen::Vector3d{0.0, 0.0, 1.0});

  mesh::Tetrahedron &t0 = sendMesh.createTetrahedron(v0, v1, v2, v3);

  // Create mesh communicator
  auto &comm = *precice::utils::IntraComm::getCommunication();

  if (context.isPrimary()) {
    com::broadcastSendMesh(comm, sendMesh);
  } else {
    mesh::Mesh recvMesh("Received Mesh", dim, testing::nextMeshID());
    // receiveMesh can also deal with delta meshes
    recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
    com::broadcastReceiveMesh(comm, recvMesh);
    BOOST_TEST(recvMesh.nVertices() == 5); // 4 + 1
    BOOST_TEST(testing::equals(recvMesh.vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 9)));
    BOOST_TEST(recvMesh.tetrahedra().size() == 1);
    BOOST_TEST(testing::equals(recvMesh.tetrahedra()[0].vertex(0).getCoords(), Eigen::Vector3d{0.0, 0.0, 0.0}));
    BOOST_TEST(recvMesh.tetrahedra()[0] == t0);
  }
}

BOOST_AUTO_TEST_SUITE_END() // Mesh
BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI
