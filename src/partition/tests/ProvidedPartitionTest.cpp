#ifndef PRECICE_NO_MPI
#include "testing/Fixtures.hpp"
#include "testing/Testing.hpp"

#include "com/CommunicateBoundingBox.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/SocketCommunication.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/M2N.hpp"
#include "m2n/PointToPointComFactory.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "partition/ProvidedPartition.hpp"
#include "partition/ReceivedPartition.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace partition;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ProvidedPartitionTests)

void setupParallelEnvironment(m2n::PtrM2N m2n)
{
  BOOST_TEST(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication masterSlaveCom = com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication   = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0) { //NASTIN
    utils::Parallel::splitCommunicator("Fluid");
    m2n->acceptMasterConnection("Fluid", "SolidMaster");
  } else if (utils::Parallel::getProcessRank() == 1) { //Master
    utils::Parallel::splitCommunicator("SolidMaster");
    m2n->requestMasterConnection("Fluid", "SolidMaster");
    utils::MasterSlave::configure(0, 3);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
    utils::Parallel::splitCommunicator("SolidSlaves");
    utils::MasterSlave::configure(1, 3);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
    utils::Parallel::splitCommunicator("SolidSlaves");
    utils::MasterSlave::configure(2, 3);
  }

  if (utils::Parallel::getProcessRank() == 1) { //Master
    masterSlaveCom->acceptConnection("SolidMaster", "SolidSlaves", "Test", utils::Parallel::getProcessRank());
    masterSlaveCom->setRankOffset(1);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
    masterSlaveCom->requestConnection("SolidMaster", "SolidSlaves", "Test", 0, 2);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
    masterSlaveCom->requestConnection("SolidMaster", "SolidSlaves", "Test", 1, 2);
  }
}

// create a communciator with two participants: each one has one master and one slave ranks
// only master-master and master-slave channels are created here
// Point to Point slave-slave channels must be created in the test.
void setupM2NBaseEnvironment(m2n::PtrM2N m2n)
{
  PRECICE_ASSERT(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication masterSlaveCom = com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication   = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0) { //Master Fluid
    utils::Parallel::splitCommunicator("FluidMaster");
    m2n->acceptMasterConnection("FluidMaster", "SolidMaster");
    utils::MasterSlave::configure(0, 2);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    utils::Parallel::splitCommunicator("FluidSlave");
    utils::MasterSlave::configure(1, 2);
  } else if (utils::Parallel::getProcessRank() == 2) { //Master Solid
    utils::Parallel::splitCommunicator("SolidMaster");
    m2n->requestMasterConnection("FluidMaster", "SolidMaster");
    utils::MasterSlave::configure(0, 2);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
    utils::Parallel::splitCommunicator("SolidSlave");
    utils::MasterSlave::configure(1, 2);
  }

  if (utils::Parallel::getProcessRank() == 0) { //Master Fluid
    utils::MasterSlave::_communication->acceptConnection("FluidMaster", "FluidSlave", "Test", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave Fluid
    utils::MasterSlave::_communication->requestConnection("FluidMaster", "FluidSlave", "Test", 0, 1);
  } else if (utils::Parallel::getProcessRank() == 2) { //Master Solid
    utils::MasterSlave::_communication->acceptConnection("SolidMaster", "SolidSlave", "Test", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave Solis
    utils::MasterSlave::_communication->requestConnection("SolidMaster", "SolidSlave", "Test", 0, 1);
  }
}

void tearDownParallelEnvironment()
{
  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
}

BOOST_AUTO_TEST_CASE(TestGatherAndCommunicate2D, *testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int  dimensions  = 2;
  bool flipNormals = false;

  if (utils::Parallel::getProcessRank() == 0) { //NASTIN
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    double safetyFactor = 0.1;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_SLAVES, safetyFactor);
    part.addM2N(m2n);
    part.communicate();

    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 4);

    for (int i = 0; i < 6; i++) {
      BOOST_TEST(pSolidzMesh->vertices()[i].getGlobalIndex() == i);
    }
  } else { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (utils::Parallel::getProcessRank() == 1) { //Master
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0;
      mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
      position << 0.0, 1.5;
      mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v1, v2);
    } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      Eigen::VectorXd position(dimensions);
      position << 0.0, 3.5;
      mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5;
      mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
      position << 0.0, 5.5;
      mesh::Vertex &v5 = pSolidzMesh->createVertex(position);
      position << 0.0, 7.0;
      mesh::Vertex &v6 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3, v4);
      pSolidzMesh->createEdge(v4, v5);
      pSolidzMesh->createEdge(v5, v6);
    }

    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    if (utils::Parallel::getProcessRank() == 1) { //master
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2] == 6);
    } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2] == 6);
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2] == 6);
    }
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestGatherAndCommunicate3D, *testing::OnSize(4))
{
  com::PtrCommunication                     participantCom = com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory   = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int  dimensions  = 3;
  bool flipNormals = false;

  if (utils::Parallel::getProcessRank() == 0) { //NASTIN
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    double safetyFactor = 0.1;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_SLAVES, safetyFactor);
    part.addM2N(m2n);
    part.communicate();

    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 6);
    BOOST_TEST(pSolidzMesh->triangles().size() == 2);

    for (int i = 0; i < 6; i++) {
      BOOST_TEST(pSolidzMesh->vertices()[i].getGlobalIndex() == i);
    }
  } else { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (utils::Parallel::getProcessRank() == 1) { //Master
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0, 0.0;
      mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
      position << 0.0, 1.5, 1.0;
      mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v1, v2);
    } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      Eigen::VectorXd position(dimensions);
      position << 0.0, 3.5, 0.1;
      mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5, 0.2;
      mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
      position << 0.0, 5.5, 0.8;
      mesh::Vertex &v5 = pSolidzMesh->createVertex(position);
      position << 0.0, 7.0, 0.4;
      mesh::Vertex &v6 = pSolidzMesh->createVertex(position);
      mesh::Edge &  e1 = pSolidzMesh->createEdge(v3, v4);
      mesh::Edge &  e2 = pSolidzMesh->createEdge(v4, v5);
      mesh::Edge &  e3 = pSolidzMesh->createEdge(v5, v3);
      mesh::Edge &  e4 = pSolidzMesh->createEdge(v3, v6);
      mesh::Edge &  e5 = pSolidzMesh->createEdge(v6, v5);

      pSolidzMesh->createTriangle(e1, e2, e3);
      pSolidzMesh->createTriangle(e4, e5, e3);
    }

    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    if (utils::Parallel::getProcessRank() == 1) { //master
      BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 3);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2] == 6);
      BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
      BOOST_TEST(pSolidzMesh->vertices()[0].getGlobalIndex() == 0);
      BOOST_TEST(pSolidzMesh->vertices()[1].getGlobalIndex() == 1);
      BOOST_TEST(pSolidzMesh->vertices()[0].isOwner() == true);
      BOOST_TEST(pSolidzMesh->vertices()[1].isOwner() == true);
      BOOST_TEST(pSolidzMesh->getVertexDistribution()[0].size() == 2);
      BOOST_TEST(pSolidzMesh->getVertexDistribution()[1].size() == 0);
      BOOST_TEST(pSolidzMesh->getVertexDistribution()[2].size() == 4);
      BOOST_TEST(pSolidzMesh->getVertexDistribution()[0][0] == 0);
      BOOST_TEST(pSolidzMesh->getVertexDistribution()[0][1] == 1);
      BOOST_TEST(pSolidzMesh->getVertexDistribution()[2][0] == 2);
      BOOST_TEST(pSolidzMesh->getVertexDistribution()[2][1] == 3);
      BOOST_TEST(pSolidzMesh->getVertexDistribution()[2][2] == 4);
      BOOST_TEST(pSolidzMesh->getVertexDistribution()[2][3] == 5);
    } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
      BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 3);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2] == 6);
      BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 3);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[0] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[1] == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets()[2] == 6);
      BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
      BOOST_TEST(pSolidzMesh->vertices()[0].getGlobalIndex() == 2);
      BOOST_TEST(pSolidzMesh->vertices()[1].getGlobalIndex() == 3);
      BOOST_TEST(pSolidzMesh->vertices()[2].getGlobalIndex() == 4);
      BOOST_TEST(pSolidzMesh->vertices()[3].getGlobalIndex() == 5);
      BOOST_TEST(pSolidzMesh->vertices()[0].isOwner() == true);
      BOOST_TEST(pSolidzMesh->vertices()[1].isOwner() == true);
      BOOST_TEST(pSolidzMesh->vertices()[2].isOwner() == true);
      BOOST_TEST(pSolidzMesh->vertices()[3].isOwner() == true);
    }
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestOnlyDistribution2D,
                     *testing::OnSize(4) * boost::unit_test::fixture<testing::MasterComFixture>())
{
  // Create mesh object
  std::string   meshName("MyMesh");
  int           dim         = 2;
  bool          flipNormals = false; // The normals of triangles, edges, vertices
  mesh::PtrMesh pMesh(new mesh::Mesh(meshName, dim, flipNormals, testing::nextMeshID()));

  if (utils::Parallel::getProcessRank() == 0) { //Master
    Eigen::VectorXd position(dim);
    position << 0.0, 0.0;
    pMesh->createVertex(position);
    position << 1.0, 0.0;
    pMesh->createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    Eigen::VectorXd position(dim);
    position << 2.0, 0.0;
    pMesh->createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
    Eigen::VectorXd position(dim);
    position << 3.0, 0.0;
    pMesh->createVertex(position);
    position << 4.0, 0.0;
    pMesh->createVertex(position);
  }

  ProvidedPartition part(pMesh);
  part.communicate();
  part.compute();

  if (utils::Parallel::getProcessRank() == 0) { //Master
    BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 5);
    BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
    BOOST_TEST(pMesh->getVertexOffsets()[0] == 2);
    BOOST_TEST(pMesh->getVertexOffsets()[1] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[2] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[3] == 5);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == true);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 5);
    BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
    BOOST_TEST(pMesh->getVertexOffsets()[0] == 2);
    BOOST_TEST(pMesh->getVertexOffsets()[1] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[2] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[3] == 5);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 2);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
    BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 5);
    BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
    BOOST_TEST(pMesh->getVertexOffsets()[0] == 2);
    BOOST_TEST(pMesh->getVertexOffsets()[1] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[2] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[3] == 5);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
    BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 5);
    BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
    BOOST_TEST(pMesh->getVertexOffsets()[0] == 2);
    BOOST_TEST(pMesh->getVertexOffsets()[1] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[2] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[3] == 5);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 3);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 4);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == true);
  }
}

BOOST_AUTO_TEST_CASE(TestCompareBoundingBoxes2D, *testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  bool        useOnlyMasterCom = false;
  bool        useTwoLevelInit  = true;
  m2n::PtrM2N m2n              = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory, useOnlyMasterCom, useTwoLevelInit));

  setupParallelEnvironment(m2n);

  int  dimensions  = 2;
  bool flipNormals = true;

  BOOST_TEST(utils::Parallel::getCommunicatorSize() == 4);

  if (utils::Parallel::getProcessRank() != 0) { //SOLIDZ

    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (utils::Parallel::getProcessRank() == 1) { //Master
      Eigen::VectorXd position(dimensions);
      position << -1.0, 0.0;
      mesh::Vertex &v0 = pSolidzMesh->createVertex(position);
      position << 1.0, 2.0;
      mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
      position << 5.0, 3.0;
      mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v0, v1);
      pSolidzMesh->createEdge(v1, v2);
    }

    else if (utils::Parallel::getProcessRank() == 2) { //Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.0, 3.5;
      mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5;
      mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3, v4);
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      Eigen::VectorXd position(dimensions);
      position << 2.5, 5.5;
      mesh::Vertex &v5 = pSolidzMesh->createVertex(position);
      position << 4.5, 7.0;
      mesh::Vertex &v6 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v5, v6);
    }

    pSolidzMesh->computeBoundingBox();
    pSolidzMesh->computeState();

    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.compareBoundingBoxes();

    if (utils::Parallel::getProcessRank() == 1) { //Master
      BOOST_TEST(pSolidzMesh->getConnectedRanks().size() == 2);
      BOOST_TEST(pSolidzMesh->getConnectedRanks()[0] == 1);
      BOOST_TEST(pSolidzMesh->getConnectedRanks()[1] == 2);
    } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
      BOOST_TEST(pSolidzMesh->getConnectedRanks().size() == 2);
      BOOST_TEST(pSolidzMesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(pSolidzMesh->getConnectedRanks()[1] == 2);
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      BOOST_TEST(pSolidzMesh->getConnectedRanks().size() == 2);
      BOOST_TEST(pSolidzMesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(pSolidzMesh->getConnectedRanks()[1] == 1);
    }

  } else { //NASTIN

    mesh::Mesh::BoundingBoxMap receivedGlobalBB;
    mesh::BoundingBox    localBB;

    // we receive other participants communicator size
    int receivedFeedbackSize = 3;
    m2n->getMasterCommunication()->receive(receivedFeedbackSize, 0);

    for (int j = 0; j < dimensions; j++) {
      localBB.setBounds(j, -1, -1);
    }
    for (int i = 0; i < receivedFeedbackSize; i++) {
      receivedGlobalBB[i] = localBB;
    }

    // we receive golbal bounding box from othe participant!
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveBoundingBoxMap(receivedGlobalBB, 0);
    // check wether we have received the correct com size
    BOOST_TEST(receivedFeedbackSize == 3);

    //check the validity of received golbal bounding box (globalBB)
    BOOST_TEST(receivedGlobalBB[0].getData(0,1) == -1);
    BOOST_TEST(receivedGlobalBB[0].getData(0,2) == 5);
    BOOST_TEST(receivedGlobalBB[0].getData(1,1) == 0);
    BOOST_TEST(receivedGlobalBB[0].getData(1,2) == 3);
    BOOST_TEST(receivedGlobalBB[1].getData(0,1) == 0);
    BOOST_TEST(receivedGlobalBB[1].getData(0,2) == 1);
    BOOST_TEST(receivedGlobalBB[1].getData(1,1) == 3.5);
    BOOST_TEST(receivedGlobalBB[1].getData(1,2) == 4.5);
    BOOST_TEST(receivedGlobalBB[2].getData(0,1) == 2.5);
    BOOST_TEST(receivedGlobalBB[2].getData(0,2) == 4.5);
    BOOST_TEST(receivedGlobalBB[2].getData(1,1) == 5.5);
    BOOST_TEST(receivedGlobalBB[2].getData(1,2) == 7.0);

    std::vector<int> connectedRanks = {0, 1, 2};
    m2n->getMasterCommunication()->send(connectedRanks, 0);

    // construct connection map
    std::map<int, std::vector<int>> connectionMap;
    connectionMap[0].push_back(1);
    connectionMap[0].push_back(2);
    connectionMap[1].push_back(0);
    connectionMap[1].push_back(2);
    connectionMap[2].push_back(0);
    connectionMap[2].push_back(1);

    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendConnectionMap(connectionMap, 0);
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestSendBoundingBoxes3D, *testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  bool        useOnlyMasterCom = false;
  bool        useTwoLevelInit  = true;
  m2n::PtrM2N m2n              = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory, useOnlyMasterCom, useTwoLevelInit));

  setupParallelEnvironment(m2n);

  int  dimensions  = 3;
  bool flipNormals = true;

  BOOST_TEST(utils::Parallel::getCommunicatorSize() == 4);

  if (utils::Parallel::getProcessRank() != 0) { //SOLIDZ

    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (utils::Parallel::getProcessRank() == 1) { //Master
      Eigen::VectorXd position(dimensions);
      position << -1.0, 0.0, -1.0;
      mesh::Vertex &v0 = pSolidzMesh->createVertex(position);
      position << 1.0, 2.0, 1.0;
      mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
      position << 5.0, 3.0, 5.0;
      mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v0, v1);
      pSolidzMesh->createEdge(v1, v2);
    }

    else if (utils::Parallel::getProcessRank() == 2) { //Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.0, 3.5, 1.0;
      mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5, 0.0;
      mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3, v4);
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      Eigen::VectorXd position(dimensions);
      position << 2.5, 5.5, 2.5;
      mesh::Vertex &v5 = pSolidzMesh->createVertex(position);
      position << 4.5, 7.0, 4.5;
      mesh::Vertex &v6 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v5, v6);
    }

    pSolidzMesh->computeBoundingBox();

    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.compareBoundingBoxes();

  } else { //NASTIN

    mesh::Mesh::BoundingBoxMap receivedGlobalBB;
    mesh::BoundingBox    localBB;

    // we receive other participants communicator size
    int remoteParComSize = 3;
    m2n->getMasterCommunication()->receive(remoteParComSize, 0);

    for (int j = 0; j < dimensions; j++) {
      localBB.setBounds(j,-1, -1);
    }
    for (int i = 0; i < remoteParComSize; i++) {
      receivedGlobalBB[i] = localBB;
    }

    // we receive golbal bounding box from othe participant!
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveBoundingBoxMap(receivedGlobalBB, 0);

    // check wether we have received the correct com size
    BOOST_TEST(remoteParComSize == 3);

    //check the validity of received golbal bounding box (globalBB)
    BOOST_TEST(receivedGlobalBB[0].getData(0,1) == -1);
    BOOST_TEST(receivedGlobalBB[0].getData(0,2) == 5);
    BOOST_TEST(receivedGlobalBB[0].getData(1,1) == 0);
    BOOST_TEST(receivedGlobalBB[0].getData(1,2) == 3);
    BOOST_TEST(receivedGlobalBB[0].getData(2,1) == -1);
    BOOST_TEST(receivedGlobalBB[0].getData(2,2) == 5);
    BOOST_TEST(receivedGlobalBB[1].getData(0,1) == 0);
    BOOST_TEST(receivedGlobalBB[1].getData(0,2) == 1);
    BOOST_TEST(receivedGlobalBB[1].getData(1,1) == 3.5);
    BOOST_TEST(receivedGlobalBB[1].getData(1,2) == 4.5);
    BOOST_TEST(receivedGlobalBB[1].getData(2,1) == 0);
    BOOST_TEST(receivedGlobalBB[1].getData(2,2) == 1);
    BOOST_TEST(receivedGlobalBB[2].getData(0,1) == 2.5);
    BOOST_TEST(receivedGlobalBB[2].getData(0,2) == 4.5);
    BOOST_TEST(receivedGlobalBB[2].getData(1,1)== 5.5);
    BOOST_TEST(receivedGlobalBB[2].getData(1,2) == 7.0);
    BOOST_TEST(receivedGlobalBB[2].getData(2,1) == 2.5);
    BOOST_TEST(receivedGlobalBB[2].getData(2,2) == 4.5);

    //send empty dummy list of connected ranks as feedback
    std::vector<int> connectedRanksList;
    m2n->getMasterCommunication()->send(connectedRanksList, 0);
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestCommunicateLocalMeshPartitions, *testing::OnSize(4))
{
  //mesh creation
  int           dimensions       = 2;
  bool          flipNormals      = true;
  double        safetyFactor     = 0.1;
  bool          useOnlyMasterCom = false;
  bool          useTwoLevelInit  = true;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));

  // create second communicator for m2n mesh and communciation map exchange
  com::PtrCommunication                     participantsCom       = com::PtrCommunication(new com::SocketCommunication());
  com::PtrCommunicationFactory              participantComFactory = com::PtrCommunicationFactory(new com::SocketCommunicationFactory);
  m2n::DistributedComFactory::SharedPointer distributionFactory   = m2n::DistributedComFactory::SharedPointer(new m2n::PointToPointComFactory(participantComFactory));
  m2n::PtrM2N                               m2n                   = m2n::PtrM2N(new m2n::M2N(participantsCom, distributionFactory, useOnlyMasterCom, useTwoLevelInit));

  setupM2NBaseEnvironment(m2n);

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    Eigen::VectorXd position(dimensions);
    position << 0.5, 0.0;
    mesh::Vertex &v1 = mesh->createVertex(position);
    position << 1.5, 0.0;
    mesh::Vertex &v2 = mesh->createVertex(position);
    position << 2.0, 1.0;
    mesh::Vertex &v3 = mesh->createVertex(position);
    position << 0.5, 1.0;
    mesh::Vertex &v4 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);
    mesh->createEdge(v2, v3);
    mesh->createEdge(v3, v4);
    mesh->createEdge(v4, v1);

    mesh->getConnectedRanks().push_back(0);

    break;
  }
  case 1: {
    Eigen::VectorXd position(dimensions);
    position << 2.5, 0.0;
    mesh::Vertex &v1 = mesh->createVertex(position);
    position << 3.5, 0.0;
    mesh::Vertex &v2 = mesh->createVertex(position);
    position << 3.5, 1.0;
    mesh::Vertex &v3 = mesh->createVertex(position);
    position << 2.0, 1.0;
    mesh::Vertex &v4 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);
    mesh->createEdge(v2, v3);
    mesh->createEdge(v3, v4);
    mesh->createEdge(v4, v1);

    mesh->getConnectedRanks().push_back(1);

    break;
  }
  case 2: {

    mesh->getConnectedRanks().push_back(0);

    break;
  }
  case 3: {

    mesh->getConnectedRanks().push_back(1);

    break;
  }
  }

  if (utils::Parallel::getProcessRank() < 2) {
    m2n->createDistributedCommunication(mesh);
    ProvidedPartition part(mesh);
    m2n->acceptSlavesPreConnection("SolidSlaves", "FluidSlaves");
    part.addM2N(m2n);
    part.communicate();
  } else {
    m2n->createDistributedCommunication(mesh);
    ReceivedPartition part(mesh, ReceivedPartition::ON_SLAVES, safetyFactor);
    m2n->requestSlavesPreConnection("SolidSlaves", "FluidSlaves");
    part.addM2N(m2n);

    part.communicate();

    BOOST_TEST(mesh->vertices().size() == 4);

    if (utils::Parallel::getProcessRank() == 2) {
      BOOST_TEST(mesh->vertices()[0].getCoords()[0] == 0.5);
      BOOST_TEST(mesh->vertices()[0].getCoords()[1] == 0.0);
      BOOST_TEST(mesh->vertices()[1].getCoords()[0] == 1.5);
      BOOST_TEST(mesh->vertices()[1].getCoords()[1] == 0.0);
      BOOST_TEST(mesh->vertices()[2].getCoords()[0] == 2.0);
      BOOST_TEST(mesh->vertices()[2].getCoords()[1] == 1.0);
      BOOST_TEST(mesh->vertices()[3].getCoords()[0] == 0.5);
      BOOST_TEST(mesh->vertices()[3].getCoords()[1] == 1.0);
    } else {
      BOOST_TEST(mesh->vertices()[0].getCoords()[0] == 2.5);
      BOOST_TEST(mesh->vertices()[0].getCoords()[1] == 0.0);
      BOOST_TEST(mesh->vertices()[1].getCoords()[0] == 3.5);
      BOOST_TEST(mesh->vertices()[1].getCoords()[1] == 0.0);
      BOOST_TEST(mesh->vertices()[2].getCoords()[0] == 3.5);
      BOOST_TEST(mesh->vertices()[2].getCoords()[1] == 1.0);
      BOOST_TEST(mesh->vertices()[3].getCoords()[0] == 2.0);
      BOOST_TEST(mesh->vertices()[3].getCoords()[1] == 1.0);
    }
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestTwoLevelRepartitioning2D, *testing::OnSize(4))
{
  //mesh creation
  int           dimensions   = 2;
  bool          flipNormals  = true;
  double        safetyFactor = 0;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh receivedMesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    Eigen::VectorXd position(dimensions);
    position << -2.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 1.0;
    mesh->createVertex(position);
    position << -1.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, 2.0;
    mesh->createVertex(position);
    position << -1.0, 2.0;
    mesh->createVertex(position);
    position << 0.0, 2.0;
    mesh->createVertex(position);

    break;
  }
  case 1: {
    Eigen::VectorXd position(dimensions);
    position << -0.5, 0.0;
    mesh->createVertex(position);
    position << 1.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, 2.0;
    mesh->createVertex(position);
    position << 2.0, 2.0;
    mesh->createVertex(position);

    break;
  }
  case 2: {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, -1.0;
    mesh->createVertex(position);
    position << -1.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, -1.0;
    mesh->createVertex(position);
    position << -2.0, -0.0;
    mesh->createVertex(position);
    position << -2.0, -1.0;
    mesh->createVertex(position);

    break;
  }
  case 3: {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, -1.0;
    mesh->createVertex(position);
    position << 1.0, -1.0;
    mesh->createVertex(position);
    position << 2.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, -1.0;
    mesh->createVertex(position);
    break;
  }
  }

  mesh->computeState();
  mesh->computeBoundingBox();

  // create the communicator for m2n mesh and communciation map exchange
  com::PtrCommunication                     participantsCom       = com::PtrCommunication(new com::SocketCommunication());
  com::PtrCommunicationFactory              participantComFactory = com::PtrCommunicationFactory(new com::SocketCommunicationFactory);
  m2n::DistributedComFactory::SharedPointer distributionFactory   = m2n::DistributedComFactory::SharedPointer(new m2n::PointToPointComFactory(participantComFactory));
  bool                                      useOnlyMasterCom      = false;
  bool                                      useTwoLevelInit       = true;
  m2n::PtrM2N                               m2n                   = m2n::PtrM2N(new m2n::M2N(participantsCom, distributionFactory, useOnlyMasterCom, useTwoLevelInit));

  setupM2NBaseEnvironment(m2n);

  if (utils::Parallel::getProcessRank() < 2) {
    m2n->createDistributedCommunication(mesh);
    ProvidedPartition part(mesh);
    part.addM2N(m2n);

    part.compareBoundingBoxes();

    if (utils::Parallel::getProcessRank() == 0) {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(mesh->getConnectedRanks()[1] == 1);
    } else {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(mesh->getConnectedRanks()[1] == 1);
    }

    m2n->acceptSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();

    if (utils::Parallel::getProcessRank() == 0) {
      BOOST_TEST(mesh->getCommunicationMap()[0][0] == 0);
      BOOST_TEST(mesh->getCommunicationMap()[0][1] == 1);
    } else {
      BOOST_TEST(mesh->getCommunicationMap()[0][0] == 0);
      BOOST_TEST(mesh->getCommunicationMap()[1][0] == 1);
      BOOST_TEST(mesh->getCommunicationMap()[1][1] == 2);
    }

  } else {
    m2n->createDistributedCommunication(receivedMesh);
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping   = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(receivedMesh, mesh);
    boundingToMapping->setMeshes(mesh, receivedMesh);

    ReceivedPartition part(receivedMesh, ReceivedPartition::ON_SLAVES, safetyFactor);

    part.addM2N(m2n);

    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);

    part.compareBoundingBoxes();

    m2n->requestSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestTwoLevelRepartitioning3D, *testing::OnSize(4))
{

  //mesh creation
  int           dimensions   = 3;
  bool          flipNormals  = true;
  double        safetyFactor = 0.0;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh receivedMesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    Eigen::VectorXd position(dimensions);
    position << -2.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, 1.0, 1.0;
    mesh->createVertex(position);
    position << -1.0, 1.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, 1.0, 1.0;
    mesh->createVertex(position);

    break;
  }
  case 1: {
    Eigen::VectorXd position(dimensions);
    position << -0.5, 0.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, 1.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, 1.0, 1.0;
    mesh->createVertex(position);
    break;
  }
  case 2: {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << -1.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << -1.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << -2.0, -0.0, 0.0;
    mesh->createVertex(position);
    position << -2.0, -1.0, 1.0;
    mesh->createVertex(position);

    break;
  }
  case 3: {
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 1.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 0.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << 1.0, -1.0, 1.0;
    mesh->createVertex(position);
    position << 2.0, 0.0, 0.0;
    mesh->createVertex(position);
    position << 2.0, -1.0, 0.0;
    mesh->createVertex(position);
    break;
  }
  }

  mesh->computeState();
  mesh->computeBoundingBox();

  // create the communicator for m2n mesh and communciation map exchange
  com::PtrCommunication                     participantsCom       = com::PtrCommunication(new com::SocketCommunication());
  com::PtrCommunicationFactory              participantComFactory = com::PtrCommunicationFactory(new com::SocketCommunicationFactory);
  m2n::DistributedComFactory::SharedPointer distributionFactory   = m2n::DistributedComFactory::SharedPointer(new m2n::PointToPointComFactory(participantComFactory));
  bool                                      useOnlyMasterCom      = false;
  bool                                      useTwoLevelInit       = true;
  m2n::PtrM2N                               m2n                   = m2n::PtrM2N(new m2n::M2N(participantsCom, distributionFactory, useOnlyMasterCom, useTwoLevelInit));

  setupM2NBaseEnvironment(m2n);

  if (utils::Parallel::getProcessRank() < 2) {
    m2n->createDistributedCommunication(mesh);
    ProvidedPartition part(mesh);
    part.addM2N(m2n);

    part.compareBoundingBoxes();

    if (utils::Parallel::getProcessRank() == 0) {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(mesh->getConnectedRanks()[1] == 1);
    } else {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(mesh->getConnectedRanks()[1] == 1);
    }

    m2n->acceptSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();

    if (utils::Parallel::getProcessRank() == 0) {
      BOOST_TEST(mesh->getCommunicationMap()[0][0] == 0);
      BOOST_TEST(mesh->getCommunicationMap()[0][1] == 1);
    } else {
      BOOST_TEST(mesh->getCommunicationMap()[0][0] == 0);
      BOOST_TEST(mesh->getCommunicationMap()[1][0] == 1);
      BOOST_TEST(mesh->getCommunicationMap()[1][1] == 2);
    }
  } else {
    m2n->createDistributedCommunication(receivedMesh);
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping   = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(receivedMesh, mesh);
    boundingToMapping->setMeshes(mesh, receivedMesh);

    ReceivedPartition part(receivedMesh, ReceivedPartition::ON_SLAVES, safetyFactor);
    part.addM2N(m2n);

    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);

    part.compareBoundingBoxes();

    m2n->requestSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
