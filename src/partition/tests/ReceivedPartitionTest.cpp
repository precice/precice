#ifndef PRECICE_NO_MPI

#include "testing/Fixtures.hpp"
#include "testing/Testing.hpp"

#include "partition/ProvidedPartition.hpp"
#include "partition/ReceivedPartition.hpp"

#include "com/MPIDirectCommunication.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/M2N.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace partition;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ReceivedPartitionTests)

void setupParallelEnvironment(m2n::PtrM2N m2n)
{
  BOOST_TEST(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication masterSlaveCom = com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication   = masterSlaveCom;

  if (utils::Parallel::getProcessRank() == 0) { //SOLIDZ
    utils::Parallel::splitCommunicator("Solid");
    m2n->acceptMasterConnection("Solid", "FluidMaster");
  } else if (utils::Parallel::getProcessRank() == 1) { //Master
    utils::Parallel::splitCommunicator("FluidMaster");
    m2n->requestMasterConnection("Solid", "FluidMaster");
    utils::MasterSlave::configure(0, 3);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
    utils::Parallel::splitCommunicator("FluidSlaves");
    utils::MasterSlave::configure(1, 3);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
    utils::Parallel::splitCommunicator("FluidSlaves");
    utils::MasterSlave::configure(2, 3);
  }

  if (utils::Parallel::getProcessRank() == 1) { //Master
    masterSlaveCom->acceptConnection("FluidMaster", "FluidSlaves", "Test", utils::Parallel::getProcessRank());
    masterSlaveCom->setRankOffset(1);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
    masterSlaveCom->requestConnection("FluidMaster", "FluidSlaves", "Test", 0, 2);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
    masterSlaveCom->requestConnection("FluidMaster", "FluidSlaves", "Test", 1, 2);
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

void createSolidzMesh2D(mesh::PtrMesh pSolidzMesh)
{
  int dimensions = 2;
  BOOST_TEST(pSolidzMesh);
  BOOST_TEST(pSolidzMesh->getDimensions() == dimensions);
  Eigen::VectorXd position(dimensions);

  position << 0.0, 0.0;
  mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
  v1.setGlobalIndex(0);
  position << 0.0, 1.95;
  mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
  v2.setGlobalIndex(1);
  position << 0.0, 2.1;
  mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
  v3.setGlobalIndex(2);
  position << 0.0, 4.5;
  mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
  v4.setGlobalIndex(3);
  position << 0.0, 5.95;
  mesh::Vertex &v5 = pSolidzMesh->createVertex(position);
  v5.setGlobalIndex(4);
  position << 0.0, 6.1;
  mesh::Vertex &v6 = pSolidzMesh->createVertex(position);
  v6.setGlobalIndex(5);
  pSolidzMesh->createEdge(v1, v2);
  pSolidzMesh->createEdge(v2, v3);
  pSolidzMesh->createEdge(v3, v4);
  pSolidzMesh->createEdge(v4, v5);
  pSolidzMesh->createEdge(v5, v6);
}

void createSolidzMesh2DSmall(mesh::PtrMesh pSolidzMesh)
{
  int dimensions = 2;
  BOOST_TEST(pSolidzMesh);
  BOOST_TEST(pSolidzMesh->getDimensions() == dimensions);
  Eigen::VectorXd position(dimensions);

  position << 0.0, 0.0;
  mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
  position << 0.0, 3.0;
  mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
  position << 0.0, 6.0;
  mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
  pSolidzMesh->createEdge(v1, v2);
  pSolidzMesh->createEdge(v2, v3);
}

void createNastinMesh2D(mesh::PtrMesh pNastinMesh)
{
  int dimensions = 2;
  BOOST_TEST(pNastinMesh);
  BOOST_TEST(pNastinMesh->getDimensions() == dimensions);

  if (utils::Parallel::getProcessRank() == 1) {

    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    pNastinMesh->createVertex(position);
    position << 0.0, 2.0;
    pNastinMesh->createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 2) {
    // not at interface
  } else if (utils::Parallel::getProcessRank() == 3) {

    Eigen::VectorXd position(dimensions);
    position << 0.0, 4.0;
    pNastinMesh->createVertex(position);
    position << 0.0, 6.0;
    pNastinMesh->createVertex(position);
  }
}

void createSolidzMesh3D(mesh::PtrMesh pSolidzMesh)
{
  int             dimensions = 3;
  Eigen::VectorXd position(dimensions);
  BOOST_TEST(pSolidzMesh);
  BOOST_TEST(pSolidzMesh->getDimensions() == dimensions);

  position << 0.0, 0.0, -0.1;
  mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
  v1.setGlobalIndex(0);
  position << -1.0, 0.0, 0.0;
  mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
  v2.setGlobalIndex(1);
  position << 1.0, 0.0, 0.0;
  mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
  v3.setGlobalIndex(2);
  position << 0.0, -1.0, 0.0;
  mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
  v4.setGlobalIndex(3);
  position << 0.0, 1.0, 0.0;
  mesh::Vertex &v5 = pSolidzMesh->createVertex(position);
  v5.setGlobalIndex(4);
  mesh::Edge &e1 = pSolidzMesh->createEdge(v1, v2);
  mesh::Edge &e2 = pSolidzMesh->createEdge(v2, v4);
  mesh::Edge &e3 = pSolidzMesh->createEdge(v4, v1);
  mesh::Edge &e4 = pSolidzMesh->createEdge(v1, v3);
  mesh::Edge &e5 = pSolidzMesh->createEdge(v3, v5);
  mesh::Edge &e6 = pSolidzMesh->createEdge(v5, v1);
  pSolidzMesh->createTriangle(e1, e2, e3);
  pSolidzMesh->createTriangle(e4, e5, e6);
}

void createNastinMesh3D(mesh::PtrMesh pNastinMesh)
{
  int dimensions = 3;
  BOOST_TEST(pNastinMesh);
  BOOST_TEST(pNastinMesh->getDimensions() == dimensions);

  if (utils::Parallel::getProcessRank() == 1) { //Master

    Eigen::VectorXd position(dimensions);
    position << -1.0, -1.0, 0.0;
    pNastinMesh->createVertex(position);
    position << -0.75, -0.75, 0.5;
    pNastinMesh->createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
    // slave1 not at interface
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave2

    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, -1.0;
    pNastinMesh->createVertex(position);
    position << 0.5, 0.5, 0.0;
    pNastinMesh->createVertex(position);
  }
}

BOOST_AUTO_TEST_CASE(RePartitionNNBroadcastFilter2D, *testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int             dimensions  = 2;
  bool            flipNormals = false;
  Eigen::VectorXd offset      = Eigen::VectorXd::Zero(dimensions);

  if (utils::Parallel::getProcessRank() == 0) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh);
    pNastinMesh->computeState();

    double safetyFactor = 0.1;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::FILTER_FIRST, safetyFactor);
    part.addM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if (utils::Parallel::getProcessRank() == 1) { //Master
      BOOST_TEST(pSolidzMesh->vertices().size() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
    } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
      BOOST_TEST(pSolidzMesh->vertices().size() == 0);
      BOOST_TEST(pSolidzMesh->edges().size() == 0);
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      BOOST_TEST(pSolidzMesh->vertices().size() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
    }
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(RePartitionNNDoubleNode2D, *testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int             dimensions  = 2;
  bool            flipNormals = false;
  Eigen::VectorXd offset      = Eigen::VectorXd::Zero(dimensions);

  if (utils::Parallel::getProcessRank() == 0) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2DSmall(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh);
    pNastinMesh->computeState();

    double safetyFactor = 0.5;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::BROADCAST_FILTER, safetyFactor);
    part.addM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if (utils::Parallel::getProcessRank() == 1) { //Master
      BOOST_TEST(pSolidzMesh->vertices().size() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
    } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
      BOOST_TEST(pSolidzMesh->vertices().size() == 0);
      BOOST_TEST(pSolidzMesh->edges().size() == 0);
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      BOOST_TEST(pSolidzMesh->vertices().size() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
    }
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(RePartitionNPPreFilterPostFilter2D, *testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int  dimensions  = 2;
  bool flipNormals = false;

  if (utils::Parallel::getProcessRank() == 0) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestProjectionMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh);

    pNastinMesh->computeState();
    double            safetyFactor = 0.1;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::FILTER_FIRST, safetyFactor);
    part.addM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if (utils::Parallel::getProcessRank() == 1) { //Master
      BOOST_TEST(pSolidzMesh->vertices().size() == 3);
      BOOST_TEST(pSolidzMesh->edges().size() == 2);
    } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
      BOOST_TEST(pSolidzMesh->vertices().size() == 0);
      BOOST_TEST(pSolidzMesh->edges().size() == 0);
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      BOOST_TEST(pSolidzMesh->vertices().size() == 3);
      BOOST_TEST(pSolidzMesh->edges().size() == 2);
    }
  }
  tearDownParallelEnvironment();
}

#ifndef PRECICE_NO_PETSC
BOOST_AUTO_TEST_CASE(RePartitionRBFGlobal2D,
                     *testing::OnSize(4) * boost::unit_test::fixture<testing::MasterComFixture>() * testing::Deleted())
{
  int           dimensions  = 2;
  bool          flipNormals = false;
  mesh::PtrMesh pMesh(new mesh::Mesh("MyMesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh pOtherMesh(new mesh::Mesh("OtherMesh", dimensions, flipNormals, testing::nextMeshID()));

  mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
      new mapping::PetRadialBasisFctMapping<mapping::ThinPlateSplines>(mapping::Mapping::CONSISTENT, dimensions,
                                                                       mapping::ThinPlateSplines(), false, false, false));
  mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
      new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
  boundingFromMapping->setMeshes(pMesh, pOtherMesh);
  boundingToMapping->setMeshes(pOtherMesh, pMesh);

  if (utils::Parallel::getProcessRank() == 0) { //Master
    createSolidzMesh2D(pMesh);
  } else { //Slaves
    createNastinMesh2D(pOtherMesh);
  }

  pOtherMesh->computeState();
  double            safetyFactor = 20.0;
  ReceivedPartition part(pMesh, ReceivedPartition::NO_FILTER, safetyFactor);
  part.setFromMapping(boundingFromMapping);
  part.setToMapping(boundingToMapping);
  part.compute();

  BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
  BOOST_TEST(pMesh->getVertexOffsets()[0] == 0);
  BOOST_TEST(pMesh->getVertexOffsets()[1] == 6);
  BOOST_TEST(pMesh->getVertexOffsets()[2] == 6);
  BOOST_TEST(pMesh->getVertexOffsets()[3] == 12);
  BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 6);

  // check if the sending and filtering worked right
  if (utils::Parallel::getProcessRank() == 0) { //Master
    BOOST_TEST(pMesh->vertices().size() == 0);
    BOOST_TEST(pMesh->edges().size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[0].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[1].size() == 6);
    BOOST_TEST(pMesh->getVertexDistribution()[2].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[3].size() == 6);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    BOOST_TEST(pMesh->vertices().size() == 6);
    BOOST_TEST(pMesh->edges().size() == 5);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[2].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[3].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[4].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[5].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[2].getGlobalIndex() == 2);
    BOOST_TEST(pMesh->vertices()[3].getGlobalIndex() == 3);
    BOOST_TEST(pMesh->vertices()[4].getGlobalIndex() == 4);
    BOOST_TEST(pMesh->vertices()[5].getGlobalIndex() == 5);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
    BOOST_TEST(pMesh->vertices().size() == 0);
    BOOST_TEST(pMesh->edges().size() == 0);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
    BOOST_TEST(pMesh->vertices().size() == 6);
    BOOST_TEST(pMesh->edges().size() == 5);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[2].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[3].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[4].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[5].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[2].getGlobalIndex() == 2);
    BOOST_TEST(pMesh->vertices()[3].getGlobalIndex() == 3);
    BOOST_TEST(pMesh->vertices()[4].getGlobalIndex() == 4);
    BOOST_TEST(pMesh->vertices()[5].getGlobalIndex() == 5);
  }
}

BOOST_AUTO_TEST_CASE(RePartitionRBFLocal2D1,
                     *testing::OnSize(4) * boost::unit_test::fixture<testing::MasterComFixture>() * testing::Deleted())
{
  int           dimensions  = 2;
  bool          flipNormals = false;
  mesh::PtrMesh pMesh(new mesh::Mesh("MyMesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh pOtherMesh(new mesh::Mesh("OtherMesh", dimensions, flipNormals, testing::nextMeshID()));

  double supportRadius = 0.25;

  mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
      new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSISTENT, dimensions,
                                                                                mapping::CompactThinPlateSplinesC2(supportRadius), false, false, false));
  mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
      new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
  boundingFromMapping->setMeshes(pMesh, pOtherMesh);
  boundingToMapping->setMeshes(pOtherMesh, pMesh);

  if (utils::Parallel::getProcessRank() == 0) { //Master
    createSolidzMesh2D(pMesh);
  } else { //Slaves
    createNastinMesh2D(pOtherMesh);
  }

  pOtherMesh->computeState();
  double            safetyFactor = 20.0;
  ReceivedPartition part(pMesh, ReceivedPartition::NO_FILTER, safetyFactor);
  part.setFromMapping(boundingFromMapping);
  part.setToMapping(boundingToMapping);
  part.compute();

  BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
  BOOST_TEST(pMesh->getVertexOffsets()[0] == 0);
  BOOST_TEST(pMesh->getVertexOffsets()[1] == 3);
  BOOST_TEST(pMesh->getVertexOffsets()[2] == 3);
  BOOST_TEST(pMesh->getVertexOffsets()[3] == 6);
  BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 6);

  // check if the sending and filtering worked right
  if (utils::Parallel::getProcessRank() == 0) { //Master
    BOOST_TEST(pMesh->vertices().size() == 0);
    BOOST_TEST(pMesh->edges().size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[0].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[1].size() == 3);
    BOOST_TEST(pMesh->getVertexDistribution()[2].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[3].size() == 3);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    BOOST_TEST(pMesh->vertices().size() == 3);
    BOOST_TEST(pMesh->edges().size() == 2);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[2].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[2].getGlobalIndex() == 2);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
    BOOST_TEST(pMesh->vertices().size() == 0);
    BOOST_TEST(pMesh->edges().size() == 0);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
    BOOST_TEST(pMesh->vertices().size() == 3);
    BOOST_TEST(pMesh->edges().size() == 2);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[2].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 3);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 4);
    BOOST_TEST(pMesh->vertices()[2].getGlobalIndex() == 5);
  }
}

BOOST_AUTO_TEST_CASE(RePartitionRBFLocal2D2,
                     *testing::OnSize(4) * boost::unit_test::fixture<testing::MasterComFixture>() * testing::Deleted())
{
  int           dimensions  = 2;
  bool          flipNormals = false;
  mesh::PtrMesh pMesh(new mesh::Mesh("MyMesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh pOtherMesh(new mesh::Mesh("OtherMesh", dimensions, flipNormals, testing::nextMeshID()));

  double supportRadius = 2.45;

  mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
      new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSISTENT, dimensions,
                                                                                mapping::CompactThinPlateSplinesC2(supportRadius), false, false, false));
  mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
      new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
  boundingFromMapping->setMeshes(pMesh, pOtherMesh);
  boundingToMapping->setMeshes(pOtherMesh, pMesh);

  if (utils::Parallel::getProcessRank() == 0) { //Master
    createSolidzMesh2D(pMesh);
  } else { //Slaves
    createNastinMesh2D(pOtherMesh);
  }

  pOtherMesh->computeState();
  double            safetyFactor = 20.0;
  ReceivedPartition part(pMesh, ReceivedPartition::NO_FILTER, safetyFactor);
  part.setFromMapping(boundingFromMapping);
  part.setToMapping(boundingToMapping);
  part.compute();

  BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
  BOOST_TEST(pMesh->getVertexOffsets()[0] == 0);
  BOOST_TEST(pMesh->getVertexOffsets()[1] == 4);
  BOOST_TEST(pMesh->getVertexOffsets()[2] == 4);
  BOOST_TEST(pMesh->getVertexOffsets()[3] == 9);
  BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 6);

  // check if the sending and filtering worked right
  if (utils::Parallel::getProcessRank() == 0) { //Master
    BOOST_TEST(pMesh->vertices().size() == 0);
    BOOST_TEST(pMesh->edges().size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[0].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[1].size() == 4);
    BOOST_TEST(pMesh->getVertexDistribution()[2].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[3].size() == 5);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    BOOST_TEST(pMesh->vertices().size() == 4);
    BOOST_TEST(pMesh->edges().size() == 3);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[2].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[3].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[2].getGlobalIndex() == 2);
    BOOST_TEST(pMesh->vertices()[3].getGlobalIndex() == 3);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
    BOOST_TEST(pMesh->vertices().size() == 0);
    BOOST_TEST(pMesh->edges().size() == 0);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
    BOOST_TEST(pMesh->vertices().size() == 5);
    BOOST_TEST(pMesh->edges().size() == 4);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[2].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[3].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[4].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 2);
    BOOST_TEST(pMesh->vertices()[2].getGlobalIndex() == 3);
    BOOST_TEST(pMesh->vertices()[3].getGlobalIndex() == 4);
    BOOST_TEST(pMesh->vertices()[4].getGlobalIndex() == 5);
  }
}

BOOST_AUTO_TEST_CASE(RePartitionRBFLocal3D,
                     *testing::OnSize(4) * boost::unit_test::fixture<testing::MasterComFixture>() * testing::Deleted())
{
  int           dimensions  = 3;
  bool          flipNormals = false;
  mesh::PtrMesh pMesh(new mesh::Mesh("MyMesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh pOtherMesh(new mesh::Mesh("OtherMesh", dimensions, flipNormals, testing::nextMeshID()));

  double supportRadius1 = 1.2;
  double supportRadius2 = 0.2;

  mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
      new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSISTENT, dimensions,
                                                                                mapping::CompactThinPlateSplinesC2(supportRadius1), false, false, false));
  mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
      new mapping::PetRadialBasisFctMapping<mapping::CompactThinPlateSplinesC2>(mapping::Mapping::CONSERVATIVE, dimensions,
                                                                                mapping::CompactThinPlateSplinesC2(supportRadius2), false, false, false));
  boundingFromMapping->setMeshes(pMesh, pOtherMesh);
  boundingToMapping->setMeshes(pOtherMesh, pMesh);

  if (utils::Parallel::getProcessRank() == 0) { //Master
    createSolidzMesh3D(pMesh);
  } else { //Slaves
    createNastinMesh3D(pOtherMesh);
  }

  pOtherMesh->computeState();
  double            safetyFactor = 20.0;
  ReceivedPartition part(pMesh, ReceivedPartition::NO_FILTER, safetyFactor);
  part.setFromMapping(boundingFromMapping);
  part.setToMapping(boundingToMapping);
  part.compute();

  BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
  BOOST_TEST(pMesh->getVertexOffsets()[0] == 0);
  BOOST_TEST(pMesh->getVertexOffsets()[1] == 5);
  BOOST_TEST(pMesh->getVertexOffsets()[2] == 5);
  BOOST_TEST(pMesh->getVertexOffsets()[3] == 10);
  BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 5);

  // check if the sending and filtering worked right
  if (utils::Parallel::getProcessRank() == 0) { //Master
    BOOST_TEST(pMesh->vertices().size() == 0);
    BOOST_TEST(pMesh->edges().size() == 0);
    BOOST_TEST(pMesh->triangles().size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[0].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[1].size() == 5);
    BOOST_TEST(pMesh->getVertexDistribution()[2].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[3].size() == 5);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    BOOST_TEST(pMesh->vertices().size() == 5);
    BOOST_TEST(pMesh->edges().size() == 6);
    BOOST_TEST(pMesh->triangles().size() == 2);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[2].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[3].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[4].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[2].getGlobalIndex() == 2);
    BOOST_TEST(pMesh->vertices()[3].getGlobalIndex() == 3);
    BOOST_TEST(pMesh->vertices()[4].getGlobalIndex() == 4);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
    BOOST_TEST(pMesh->vertices().size() == 0);
    BOOST_TEST(pMesh->edges().size() == 0);
    BOOST_TEST(pMesh->triangles().size() == 0);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
    BOOST_TEST(pMesh->vertices().size() == 5);
    BOOST_TEST(pMesh->edges().size() == 6);
    BOOST_TEST(pMesh->triangles().size() == 2);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == false);
    BOOST_TEST(pMesh->vertices()[2].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[3].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[4].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[2].getGlobalIndex() == 2);
    BOOST_TEST(pMesh->vertices()[3].getGlobalIndex() == 3);
    BOOST_TEST(pMesh->vertices()[4].getGlobalIndex() == 4);
  }
}

#endif // PRECICE_NO_PETSC

BOOST_AUTO_TEST_CASE(RePartitionNPBroadcastFilter3D, *testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int  dimensions  = 3;
  bool flipNormals = false;

  if (utils::Parallel::getProcessRank() == 0) { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh3D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
  } else {
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestProjectionMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh3D(pNastinMesh);

    pNastinMesh->computeState();
    double            safetyFactor = 20.0;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::FILTER_FIRST, safetyFactor);
    part.addM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicate();
    part.compute();

    // check if the sending and filtering worked right
    if (utils::Parallel::getProcessRank() == 1) { //Master
      BOOST_TEST(pSolidzMesh->vertices().size() == 2);
      BOOST_TEST(pSolidzMesh->edges().size() == 1);
      BOOST_TEST(pSolidzMesh->triangles().size() == 0);
    } else if (utils::Parallel::getProcessRank() == 2) { //Slave1
      BOOST_TEST(pSolidzMesh->vertices().size() == 0);
      BOOST_TEST(pSolidzMesh->edges().size() == 0);
      BOOST_TEST(pSolidzMesh->triangles().size() == 0);
    } else if (utils::Parallel::getProcessRank() == 3) { //Slave2
      BOOST_TEST(pSolidzMesh->vertices().size() == 3);
      BOOST_TEST(pSolidzMesh->edges().size() == 3);
      BOOST_TEST(pSolidzMesh->triangles().size() == 1);
    }
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestRepartitionAndDistribution2D,
                     *testing::OnSize(4) * boost::unit_test::fixture<testing::MasterComFixture>())
{
  // Create mesh object
  int           dimensions  = 2;
  bool          flipNormals = false; // The normals of triangles, edges, vertices
  mesh::PtrMesh pMesh(new mesh::Mesh("MyMesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh pOtherMesh(new mesh::Mesh("OtherMesh", dimensions, flipNormals, testing::nextMeshID()));

  mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
      new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
  boundingFromMapping->setMeshes(pMesh, pOtherMesh);

  if (utils::Parallel::getProcessRank() == 0) { //Master
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    mesh::Vertex &v1 = pMesh->createVertex(position);
    v1.setGlobalIndex(0);
    position << 1.0, 0.0;
    mesh::Vertex &v2 = pMesh->createVertex(position);
    v2.setGlobalIndex(1);
    position << 2.0, 0.0;
    mesh::Vertex &v3 = pMesh->createVertex(position);
    v3.setGlobalIndex(2);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    pOtherMesh->createVertex(position);
    position << 0.8, 0.0;
    pOtherMesh->createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
    Eigen::VectorXd position(dimensions);
    position << 1.0, 0.0;
    pOtherMesh->createVertex(position);
    position << 1.2, 0.0;
    pOtherMesh->createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
  }

  pOtherMesh->computeState();
  double            safetyFactor = 20.0; //should not filter out anything here
  ReceivedPartition part(pMesh, ReceivedPartition::FILTER_FIRST, safetyFactor);
  part.setFromMapping(boundingFromMapping);
  part.compute();

  if (utils::Parallel::getProcessRank() == 0) { //Master
    BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 3);
    BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
    BOOST_TEST(pMesh->getVertexOffsets()[0] == 0);
    BOOST_TEST(pMesh->getVertexOffsets()[1] == 2);
    BOOST_TEST(pMesh->getVertexOffsets()[2] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[3] == 3);
    BOOST_TEST(pMesh->getVertexDistribution()[0].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[1].size() == 2);
    BOOST_TEST(pMesh->getVertexDistribution()[2].size() == 1);
    BOOST_TEST(pMesh->getVertexDistribution()[3].size() == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[1][0] == 0);
    BOOST_TEST(pMesh->getVertexDistribution()[1][1] == 1);
    BOOST_TEST(pMesh->getVertexDistribution()[2][0] == 1);
    BOOST_TEST(pMesh->vertices().size() == 0);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 3);
    BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
    BOOST_TEST(pMesh->getVertexOffsets()[0] == 0);
    BOOST_TEST(pMesh->getVertexOffsets()[1] == 2);
    BOOST_TEST(pMesh->getVertexOffsets()[2] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[3] == 3);
    BOOST_TEST(pMesh->vertices().size() == 2);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pMesh->vertices()[1].isOwner() == false);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
    BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 3);
    BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
    BOOST_TEST(pMesh->getVertexOffsets()[0] == 0);
    BOOST_TEST(pMesh->getVertexOffsets()[1] == 2);
    BOOST_TEST(pMesh->getVertexOffsets()[2] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[3] == 3);
    BOOST_TEST(pMesh->vertices().size() == 1);
    BOOST_TEST(pMesh->vertices()[0].getGlobalIndex() == 1);
    BOOST_TEST(pMesh->vertices()[0].isOwner() == true);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
    BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 3);
    BOOST_TEST(pMesh->getVertexOffsets().size() == 4);
    BOOST_TEST(pMesh->getVertexOffsets()[0] == 0);
    BOOST_TEST(pMesh->getVertexOffsets()[1] == 2);
    BOOST_TEST(pMesh->getVertexOffsets()[2] == 3);
    BOOST_TEST(pMesh->getVertexOffsets()[3] == 3);
    BOOST_TEST(pMesh->vertices().size() == 0);
  }
}

BOOST_FIXTURE_TEST_CASE(ProvideAndReceiveCouplingMode, testing::M2NFixture,
                        *testing::MinRanks(2) * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2)
    return;

  int  dimensions  = 2;
  bool flipNormals = false;

  if (utils::Parallel::getProcessRank() == 0) {
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    createSolidzMesh2D(pSolidzMesh);
    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 5);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 1);
    BOOST_TEST(pSolidzMesh->getVertexOffsets()[0] == 6);
    BOOST_TEST(pSolidzMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pSolidzMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pSolidzMesh->vertices()[2].getGlobalIndex() == 2);
    BOOST_TEST(pSolidzMesh->vertices()[3].getGlobalIndex() == 3);
    BOOST_TEST(pSolidzMesh->vertices()[4].getGlobalIndex() == 4);
    BOOST_TEST(pSolidzMesh->vertices()[5].getGlobalIndex() == 5);
    BOOST_TEST(pSolidzMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[1].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[2].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[3].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[4].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[5].isOwner() == true);
  } else {
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mesh::PtrMesh       pOtherMesh(new mesh::Mesh("OtherMesh", dimensions, flipNormals, testing::nextMeshID()));
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh, pOtherMesh);

    double            safetyFactor = 0.1;
    ReceivedPartition part(pSolidzMesh, ReceivedPartition::FILTER_FIRST, safetyFactor);
    part.setFromMapping(boundingFromMapping);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 5);
    BOOST_TEST(pSolidzMesh->getVertexOffsets().size() == 1);
    BOOST_TEST(pSolidzMesh->getVertexOffsets()[0] == 6);
    BOOST_TEST(pSolidzMesh->vertices()[0].getGlobalIndex() == 0);
    BOOST_TEST(pSolidzMesh->vertices()[1].getGlobalIndex() == 1);
    BOOST_TEST(pSolidzMesh->vertices()[2].getGlobalIndex() == 2);
    BOOST_TEST(pSolidzMesh->vertices()[3].getGlobalIndex() == 3);
    BOOST_TEST(pSolidzMesh->vertices()[4].getGlobalIndex() == 4);
    BOOST_TEST(pSolidzMesh->vertices()[5].getGlobalIndex() == 5);
    BOOST_TEST(pSolidzMesh->vertices()[0].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[1].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[2].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[3].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[4].isOwner() == true);
    BOOST_TEST(pSolidzMesh->vertices()[5].isOwner() == true);
  }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
