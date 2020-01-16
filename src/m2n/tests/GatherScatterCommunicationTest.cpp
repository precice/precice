#ifndef PRECICE_NO_MPI

#include "com/MPIDirectCommunication.hpp"
#include "m2n/DistributedComFactory.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/Testing.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

BOOST_AUTO_TEST_SUITE(M2NTests)

using namespace precice;
using namespace m2n;

BOOST_AUTO_TEST_CASE(GatherScatterTest, *testing::OnSize(4))
{
  BOOST_TEST(utils::Parallel::getCommunicatorSize() == 4);
  utils::MasterSlave::reset();

  com::PtrCommunication                     participantCom = com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory =
      m2n::DistributedComFactory::SharedPointer(
          new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N           m2n            = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));
  com::PtrCommunication masterSlaveCom = com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication   = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0) { // Participant 1
    utils::Parallel::splitCommunicator("Part1");
    utils::MasterSlave::configure(0, 1);
  } else if (utils::Parallel::getProcessRank() == 1) { // Participant 2 - Master
    utils::Parallel::splitCommunicator("Part2Master");
    utils::MasterSlave::configure(0, 3);
    masterSlaveCom->acceptConnection("Part2Master", "Part2Slaves", "Test", utils::Parallel::getProcessRank());
    masterSlaveCom->setRankOffset(1);
  } else if (utils::Parallel::getProcessRank() == 2) { // Participant 2 - Slave1
    utils::Parallel::splitCommunicator("Part2Slaves");
    utils::MasterSlave::configure(1, 3);
    masterSlaveCom->requestConnection("Part2Master", "Part2Slaves", "Test", 0, 2);
  } else if (utils::Parallel::getProcessRank() == 3) { // Participant 2 - Slave2
    utils::Parallel::splitCommunicator("Part2Slaves");
    utils::MasterSlave::configure(2, 3);
    masterSlaveCom->requestConnection("Part2Master", "Part2Slaves", "Test", 1, 2);
  }

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0) { // Part1
    m2n->acceptMasterConnection("Part1", "Part2Master");
  } else if (utils::Parallel::getProcessRank() == 1) { // Part2 Master
    m2n->requestMasterConnection("Part1", "Part2Master");
  } else if (utils::Parallel::getProcessRank() == 2) { // Part2 Slave1
    m2n->requestMasterConnection("Part1", "Part2Master");
  } else if (utils::Parallel::getProcessRank() == 3) { // Part2 Slave2
    m2n->requestMasterConnection("Part1", "Part2Master");
  }

  utils::Parallel::synchronizeProcesses();

  int             dimensions       = 2;
  int             numberOfVertices = 6;
  bool            flipNormals      = false;
  int             valueDimension   = 1;
  Eigen::VectorXd offset           = Eigen::VectorXd::Zero(dimensions);

  if (utils::Parallel::getProcessRank() == 0) { // Part1
    mesh::PtrMesh pMesh(new mesh::Mesh("Mesh", dimensions, flipNormals, testing::nextMeshID()));
    m2n->createDistributedCommunication(pMesh);

    pMesh->setGlobalNumberOfVertices(numberOfVertices);
    pMesh->getVertexDistribution()[0].push_back(0);
    pMesh->getVertexDistribution()[0].push_back(1);
    pMesh->getVertexDistribution()[0].push_back(2);
    pMesh->getVertexDistribution()[0].push_back(3);
    pMesh->getVertexDistribution()[0].push_back(4);
    pMesh->getVertexDistribution()[0].push_back(5);

    m2n->acceptSlavesConnection("Part1", "Part2Master");
    Eigen::VectorXd values = Eigen::VectorXd::Zero(numberOfVertices);
    values << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    m2n->send(values.data(), numberOfVertices, pMesh->getID(), valueDimension);
    m2n->receive(values.data(), numberOfVertices, pMesh->getID(),
                 valueDimension);
    BOOST_TEST(values[0] == 2.0);
    BOOST_TEST(values[1] == 4.0);
    BOOST_TEST(values[2] == 6.0);
    BOOST_TEST(values[3] == 16.0);
    BOOST_TEST(values[4] == 10.0);
    BOOST_TEST(values[5] == 12.0);

  } else {
    mesh::PtrMesh pMesh(new mesh::Mesh("Mesh", dimensions, flipNormals, testing::nextMeshID()));
    m2n->createDistributedCommunication(pMesh);
    m2n->requestSlavesConnection("Part1", "Part2Master");

    if (utils::Parallel::getProcessRank() == 1) { // Master
      pMesh->setGlobalNumberOfVertices(numberOfVertices);
      pMesh->getVertexDistribution()[0].push_back(0);
      pMesh->getVertexDistribution()[0].push_back(1);
      pMesh->getVertexDistribution()[0].push_back(3);
      pMesh->getVertexDistribution()[2].push_back(2);
      pMesh->getVertexDistribution()[2].push_back(3);
      pMesh->getVertexDistribution()[2].push_back(4);
      pMesh->getVertexDistribution()[2].push_back(5);

      Eigen::Vector3d values(0.0, 0.0, 0.0);
      m2n->receive(values.data(), 3, pMesh->getID(), valueDimension);
      BOOST_TEST(values[0] == 1.0);
      BOOST_TEST(values[1] == 2.0);
      BOOST_TEST(values[2] == 4.0);
      values = values * 2;
      m2n->send(values.data(), 3, pMesh->getID(), valueDimension);
    } else if (utils::Parallel::getProcessRank() == 2) { // Slave1
      Eigen::VectorXd values;
      m2n->receive(values.data(), 0, pMesh->getID(), valueDimension);
      m2n->send(values.data(), 0, pMesh->getID(), valueDimension);
    } else if (utils::Parallel::getProcessRank() == 3) { // Slave2
      Eigen::Vector4d values(0.0, 0.0, 0.0, 0.0);
      m2n->receive(values.data(), 4, pMesh->getID(), valueDimension);
      BOOST_TEST(values[0] == 3.0);
      BOOST_TEST(values[1] == 4.0);
      BOOST_TEST(values[2] == 5.0);
      BOOST_TEST(values[3] == 6.0);
      values = values * 2;
      m2n->send(values.data(), 4, pMesh->getID(), valueDimension);
    }
  }

  utils::MasterSlave::_communication.reset();
  utils::MasterSlave::configure(utils::Parallel::getProcessRank(), utils::Parallel::getCommunicatorSize());
  utils::MasterSlave::reset();

  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
}

BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
