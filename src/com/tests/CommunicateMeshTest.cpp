#ifndef PRECICE_NO_MPI

#include "com/CommunicateMesh.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Vertex.hpp"
#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace precice::com;

BOOST_AUTO_TEST_SUITE(CommunicationTests)

BOOST_AUTO_TEST_SUITE(MeshTests)

BOOST_AUTO_TEST_CASE(twoSolvers,
                     * testing::MinRanks(2))
{
  utils::Parallel::synchronizeProcesses();
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::PropertyContainer::resetPropertyIDCounter();

  std::string participant0("rank0");
  std::string participant1("rank1");

  for (int dim = 2; dim <= 3; dim++) {
    // Build mesh to communicate for rank0
    mesh::Mesh mesh("MyMesh", dim, false);
    if (utils::Parallel::getProcessRank() == 0) {
      mesh::Vertex &v0 = mesh.createVertex(Eigen::VectorXd::Constant(dim, 0));
      mesh::Vertex &v1 = mesh.createVertex(Eigen::VectorXd::Constant(dim, 1));
      mesh::Vertex &v2 = mesh.createVertex(Eigen::VectorXd::Constant(dim, 2));

      mesh.createEdge(v0, v1);
      mesh.createEdge(v1, v2);
      mesh.createEdge(v2, v0);
    }

    // Create mesh communicator
    std::vector<int> involvedRanks = {0, 1};
    MPI_Comm         comm          = utils::Parallel::getRestrictedCommunicator(involvedRanks);
    if (utils::Parallel::getProcessRank() < 2) {
      utils::Parallel::setGlobalCommunicator(comm);
      BOOST_TEST(utils::Parallel::getCommunicatorSize() == 2);
      com::PtrCommunication com(new com::MPIDirectCommunication());
      CommunicateMesh       comMesh(com);

      if (utils::Parallel::getProcessRank() == 0) {
        utils::Parallel::splitCommunicator(participant0);
        com->acceptConnection(participant0, participant1, utils::Parallel::getProcessRank());
        comMesh.sendMesh(mesh, 0);
        BOOST_TEST(mesh.vertices().size() == 3);
        BOOST_TEST(mesh.edges().size() == 3);
        BOOST_TEST(testing::equals(mesh.vertices()[0].getCoords(), Eigen::VectorXd::Constant(dim, 0)));
        BOOST_TEST(testing::equals(mesh.vertices()[1].getCoords(), Eigen::VectorXd::Constant(dim, 1)));
        BOOST_TEST(testing::equals(mesh.vertices()[2].getCoords(), Eigen::VectorXd::Constant(dim, 2)));
      } else if (utils::Parallel::getProcessRank() == 1) {
        mesh.createVertex(Eigen::VectorXd::Constant(dim, 9)); // new version receiveMesh can also deal with delta meshes
        utils::Parallel::splitCommunicator(participant1);
        com->requestConnection(participant0, participant1, 0, 1);
        comMesh.receiveMesh(mesh, 0);
        BOOST_TEST(mesh.vertices().size() == 4);
        BOOST_TEST(mesh.edges().size() == 3);
        BOOST_TEST(testing::equals(mesh.vertices()[0].getCoords(), Eigen::VectorXd::Constant(dim, 9)));
        BOOST_TEST(testing::equals(mesh.vertices()[1].getCoords(), Eigen::VectorXd::Constant(dim, 0)));
        BOOST_TEST(testing::equals(mesh.vertices()[2].getCoords(), Eigen::VectorXd::Constant(dim, 1)));
        BOOST_TEST(testing::equals(mesh.vertices()[3].getCoords(), Eigen::VectorXd::Constant(dim, 2)));
      }
      com->closeConnection();

      BOOST_TEST(testing::equals(mesh.edges()[0].vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 0)));
      BOOST_TEST(testing::equals(mesh.edges()[0].vertex(1).getCoords(), Eigen::VectorXd::Constant(dim, 1)));
      BOOST_TEST(testing::equals(mesh.edges()[1].vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 1)));
      BOOST_TEST(testing::equals(mesh.edges()[1].vertex(1).getCoords(), Eigen::VectorXd::Constant(dim, 2)));
      BOOST_TEST(testing::equals(mesh.edges()[2].vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 2)));
      BOOST_TEST(testing::equals(mesh.edges()[2].vertex(1).getCoords(), Eigen::VectorXd::Constant(dim, 0)));
      utils::Parallel::clearGroups();
      utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
    }
  }
}
BOOST_AUTO_TEST_SUITE_END() // Mesh

BOOST_AUTO_TEST_SUITE_END() // Communication

#endif // not PRECICE_NO_MPI
