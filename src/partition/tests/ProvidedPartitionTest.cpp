#ifndef PRECICE_NO_MPI
#include "testing/Fixtures.hpp"
#include "testing/Testing.hpp"

#include "partition/ProvidedPartition.hpp"
#include "partition/ReceivedPartition.hpp"

#include "com/MPIDirectCommunication.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/M2N.hpp"
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
    utils::MasterSlave::configure(3, 2);
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

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::BROADCAST_FILTER, safetyFactor);
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

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::BROADCAST_FILTER, safetyFactor);
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
    BOOST_TEST(pMesh->getVertexDistribution()[0].size() == 2);
    BOOST_TEST(pMesh->getVertexDistribution()[1].size() == 1);
    BOOST_TEST(pMesh->getVertexDistribution()[2].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[3].size() == 2);
    BOOST_TEST(pMesh->getVertexDistribution()[0][0] == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[0][1] == 1);
    BOOST_TEST(pMesh->getVertexDistribution()[1][0] == 2);
    BOOST_TEST(pMesh->getVertexDistribution()[3][0] == 3);
    BOOST_TEST(pMesh->getVertexDistribution()[3][1] == 4);
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

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
