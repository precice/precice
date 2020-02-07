#ifndef PRECICE_NO_MPI

#include "com/CommunicateMesh.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(VertexEdgeMesh,
                     *testing::MinRanks(2))
{
  utils::Parallel::synchronizeProcesses();
  BOOST_TEST(utils::Parallel::getCommunicatorSize() > 1);

  std::string participant0("rank0");
  std::string participant1("rank1");

  for (int dim = 2; dim <= 3; dim++) {
    mesh::Mesh    sendMesh("Sent Mesh", dim, false, testing::nextMeshID());
    mesh::Vertex &v0 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 0));
    mesh::Vertex &v1 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
    mesh::Vertex &v2 = sendMesh.createVertex(Eigen::VectorXd::Constant(dim, 2));
    mesh::Edge &  e0 = sendMesh.createEdge(v0, v1);
    mesh::Edge &  e1 = sendMesh.createEdge(v1, v2);
    mesh::Edge &  e2 = sendMesh.createEdge(v2, v0);

    // Create mesh communicator
    std::vector<int> involvedRanks = {0, 1};
    MPI_Comm         comm          = utils::Parallel::getRestrictedCommunicator(involvedRanks);

    if (utils::Parallel::getProcessRank() < 2) {
      utils::Parallel::setGlobalCommunicator(comm);
      com::PtrCommunication com(new com::MPIDirectCommunication());
      CommunicateMesh       comMesh(com);

      if (utils::Parallel::getProcessRank() == 0) {
        utils::Parallel::splitCommunicator(participant0);
        com->acceptConnection(participant0, participant1, "", utils::Parallel::getProcessRank());
        comMesh.sendMesh(sendMesh, 0);
      } else if (utils::Parallel::getProcessRank() == 1) {
        // receiveMesh can also deal with delta meshes
        mesh::Mesh recvMesh("Received Mesh", dim, false, testing::nextMeshID());
        recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
        utils::Parallel::splitCommunicator(participant1);
        com->requestConnection(participant0, participant1, "", 0, 1);
        comMesh.receiveMesh(recvMesh, 0);
        BOOST_TEST(recvMesh.vertices().size() == 4);
        BOOST_TEST(testing::equals(recvMesh.vertices()[0].getCoords(), Eigen::VectorXd::Constant(dim, 9)));
        BOOST_TEST(recvMesh.vertices()[1] == v0);
        BOOST_TEST(recvMesh.vertices()[2] == v1);
        BOOST_TEST(recvMesh.vertices()[3] == v2);
        BOOST_TEST(recvMesh.edges()[0] == e0);
        BOOST_TEST(recvMesh.edges()[1] == e1);
        BOOST_TEST(recvMesh.edges()[2] == e2);
      }
      com->closeConnection();

      utils::Parallel::clearGroups();
      utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
    }
  }
}

BOOST_AUTO_TEST_CASE(VertexEdgeTriangleMesh,
                     *testing::MinRanks(2))
{
  utils::Parallel::synchronizeProcesses();
  BOOST_TEST(utils::Parallel::getCommunicatorSize() > 1);

  std::string participant0("rank0");
  std::string participant1("rank1");

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
  std::vector<int> involvedRanks = {0, 1};
  MPI_Comm         comm          = utils::Parallel::getRestrictedCommunicator(involvedRanks);

  if (utils::Parallel::getProcessRank() < 2) {
    utils::Parallel::setGlobalCommunicator(comm);
    com::PtrCommunication com(new com::MPIDirectCommunication());
    CommunicateMesh       comMesh(com);

    if (utils::Parallel::getProcessRank() == 0) {
      utils::Parallel::splitCommunicator(participant0);
      com->acceptConnection(participant0, participant1, "", utils::Parallel::getProcessRank());
      comMesh.sendMesh(sendMesh, 0);
    } else if (utils::Parallel::getProcessRank() == 1) {
      mesh::Mesh recvMesh("Received Mesh", dim, false, testing::nextMeshID());
      // receiveMesh can also deal with delta meshes
      recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
      utils::Parallel::splitCommunicator(participant1);
      com->requestConnection(participant0, participant1, "", 0, 1);
      comMesh.receiveMesh(recvMesh, 0);
      BOOST_TEST(recvMesh.vertices().size() == 4);
      BOOST_TEST(testing::equals(recvMesh.vertices()[0].getCoords(), Eigen::VectorXd::Constant(dim, 9)));
      BOOST_TEST(recvMesh.vertices()[1] == v0);
      BOOST_TEST(recvMesh.vertices()[2] == v1);
      BOOST_TEST(recvMesh.vertices()[3] == v2);
      BOOST_TEST(recvMesh.edges()[0] == e0);
      BOOST_TEST(recvMesh.edges()[1] == e1);
      BOOST_TEST(recvMesh.edges()[2] == e2);

      BOOST_TEST(recvMesh.triangles()[0] == t0);
    }
    com->closeConnection();

    utils::Parallel::clearGroups();
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
  }
}

BOOST_AUTO_TEST_CASE(BroadcastVertexEdgeTriangleMesh,
                     *testing::MinRanks(2))
{
  utils::Parallel::synchronizeProcesses();
  BOOST_TEST(utils::Parallel::getCommunicatorSize() > 1);

  std::string participant0("rank0");
  std::string participant1("rank1");

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
  std::vector<int> involvedRanks = {0, 1};
  MPI_Comm         comm          = utils::Parallel::getRestrictedCommunicator(involvedRanks);

  if (utils::Parallel::getProcessRank() < 2) {
    utils::Parallel::setGlobalCommunicator(comm);
    com::PtrCommunication com(new com::MPIDirectCommunication());
    CommunicateMesh       comMesh(com);

    if (utils::Parallel::getProcessRank() == 0) {
      utils::Parallel::splitCommunicator(participant0);
      com->acceptConnection(participant0, participant1, "", utils::Parallel::getProcessRank());
      comMesh.broadcastSendMesh(sendMesh);
    } else if (utils::Parallel::getProcessRank() == 1) {
      mesh::Mesh recvMesh("Received Mesh", dim, false, testing::nextMeshID());
      // receiveMesh can also deal with delta meshes
      recvMesh.createVertex(Eigen::VectorXd::Constant(dim, 9));
      utils::Parallel::splitCommunicator(participant1);
      com->requestConnection(participant0, participant1, "", 0, 1);
      comMesh.broadcastReceiveMesh(recvMesh);
      BOOST_TEST(recvMesh.vertices().size() == 4);
      BOOST_TEST(testing::equals(recvMesh.vertices()[0].getCoords(), Eigen::VectorXd::Constant(dim, 9)));
      BOOST_TEST(recvMesh.vertices()[1] == v0);
      BOOST_TEST(recvMesh.vertices()[2] == v1);
      BOOST_TEST(recvMesh.vertices()[3] == v2);
      BOOST_TEST(recvMesh.edges()[0] == e0);
      BOOST_TEST(recvMesh.edges()[1] == e1);
      BOOST_TEST(recvMesh.edges()[2] == e2);
      BOOST_TEST(recvMesh.triangles()[0] == t0);
    }
    com->closeConnection();

    utils::Parallel::clearGroups();
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
  }
}

BOOST_AUTO_TEST_SUITE_END() // Mesh
BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI
