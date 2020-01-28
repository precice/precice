#ifndef PRECICE_NO_MPI
#include "testing/Fixtures.hpp"
#include "testing/Testing.hpp"

#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicationFactory.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/SocketCommunication.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "partition/ProvidedBoundingBox.hpp"
#include "partition/ReceivedBoundingBox.hpp"
#include "partition/SharedPointer.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace partition;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ReceivedBoundingBoxTests)

void setupParallelEnvironment(m2n::PtrM2N m2n, const testing::TestContext &context)
{
  // Establish and configure the master-master connection
  if (context.isNamed("Fluid")) {
    BOOST_TEST(context.hasSize(3));
    m2n->acceptMasterConnection("Fluid", "SolidMaster");
  }

  if (context.isNamed("Solid")) {
    BOOST_TEST(context.hasSize(1));
    m2n->requestMasterConnection("Fluid", "SolidMaster");
  }

  BOOST_TEST(m2n->isConnected());
}

void tearDownParallelEnvironment()
{
  mesh::Data::resetDataCount();
}

void createNastinMesh2D(mesh::PtrMesh pNastinMesh, int rank)
{
  int dimensions = 2;
  PRECICE_ASSERT(pNastinMesh.use_count() > 0);
  PRECICE_ASSERT(pNastinMesh->getDimensions() == dimensions);

  if (rank == 0) {

    Eigen::VectorXd position(dimensions);
    position << 0.10, 0.10;
    pNastinMesh->createVertex(position);
    position << 0.90, 0.90;
    pNastinMesh->createVertex(position);
  } else if (rank == 1) {
    // not at interface
  } else if (rank == 2) {

    Eigen::VectorXd position(dimensions);
    position << 2.1, 2.1;
    pNastinMesh->createVertex(position);
    position << 2.9, 2.9;
    pNastinMesh->createVertex(position);
  }
}

void createNastinMesh3D(mesh::PtrMesh pNastinMesh, int rank)
{
  int dimensions = 3;
  PRECICE_ASSERT(pNastinMesh.use_count() > 0);
  PRECICE_ASSERT(pNastinMesh->getDimensions() == dimensions);

  if (rank == 0) {

    Eigen::VectorXd position(dimensions);
    position << 0.10, 0.10, 0.1;
    pNastinMesh->createVertex(position);
    position << 0.90, 0.90, 0.9;
    pNastinMesh->createVertex(position);
  } else if (rank == 1) {
    // not at interface
  } else if (ranks == 2) {

    Eigen::VectorXd position(dimensions);
    position << 2.1, 2.1, 2.1;
    pNastinMesh->createVertex(position);
    position << 2.9, 2.9, 2.1;
    pNastinMesh->createVertex(position);
  }
}

BOOST_AUTO_TEST_CASE(TestConnectionMap2D)
{
  PRECICE_TEST("Fluid"_on(3_ranks).setupMasterSlaves(), "Solid"_on(1_rank), Require::Events);
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int  dimensions  = 2;
  bool flipNormals = true;

  // construct send global boundingbox
  mesh::Mesh::BoundingBoxMap sendGlobalBB;
  mesh::Mesh::BoundingBox    initialBB;
  for (int remoteRank = 0; remoteRank < 3; remoteRank++) {
    for (int i = 0; i < dimensions; i++) {
      initialBB.push_back(std::make_pair(3 - remoteRank - 1, 3 - remoteRank));
    }
    sendGlobalBB[remoteRank] = initialBB;
    initialBB.clear();
  }

  if (context.isNamed("Solid")) {
    std::vector<int>                connectedRanksList;
    int                             connectionMapSize = 0;
    std::map<int, std::vector<int>> receivedConnectionMap;
    mesh::PtrMesh                   pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    m2n->getMasterCommunication()->send(3, 0);
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendBoundingBoxMap(sendGlobalBB, 0);
    m2n->getMasterCommunication()->receive(connectedRanksList, 0);
    connectionMapSize = connectedRanksList.size();
    BOOST_TEST(connectionMapSize == 2);

    std::vector<int> connectedRanks;
    connectedRanks.push_back(-1);
    for (auto &rank : connectedRanksList) {
      receivedConnectionMap[rank] = connectedRanks;
    }

    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveConnectionMap(receivedConnectionMap, 0);

    // test whether we receive correct connection map
    BOOST_TEST(receivedConnectionMap[0][0] == 2);
    BOOST_TEST(receivedConnectionMap[2][0] == 0);

  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));

    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh2D(pNastinMesh);
    pNastinMesh->computeState();

    double safetyFactor = 0.0;

    ReceivedBoundingBox part(pSolidzMesh, safetyFactor);
    part.addM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicateBoundingBox();
    part.computeBoundingBox();
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestConnectionMap3D)
{
  PRECICE_TEST("Fluid"_on(3_ranks).setupMasterSlaves(), "Solid"_on(1_rank), Require::Events);
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int  dimensions  = 3;
  bool flipNormals = true;

  // construct send global boundingbox
  mesh::Mesh::BoundingBoxMap sendGlobalBB;
  mesh::Mesh::BoundingBox    initialBB;
  for (int remoteRank = 0; remoteRank < 3; remoteRank++) {
    for (int i = 0; i < dimensions; i++) {
      initialBB.push_back(std::make_pair(3 - remoteRank - 1, 3 - remoteRank));
    }
    sendGlobalBB[remoteRank] = initialBB;
    initialBB.clear();
  }

  if (context.isNamed("Solid")) {
    std::vector<int>                connectedRanksList;
    int                             connectionMapSize = 0;
    std::map<int, std::vector<int>> receivedConnectionMap;
    mesh::PtrMesh                   pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    m2n->getMasterCommunication()->send(3, 0);
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendBoundingBoxMap(sendGlobalBB, 0);
    m2n->getMasterCommunication()->receive(connectedRanksList, 0);
    connectionMapSize = connectedRanksList.size();
    BOOST_TEST(connectionMapSize == 2);

    std::vector<int> connectedRanks;
    connectedRanks.push_back(-1);
    for (auto &rank : connectedRanksList) {
      receivedConnectionMap[rank] = connectedRanks;
    }

    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveConnectionMap(receivedConnectionMap, 0);

    // test whether we receive correct connection map
    BOOST_TEST(receivedConnectionMap[0][0] == 2);
    BOOST_TEST(receivedConnectionMap[2][0] == 0);

  } else {
    BOOST_TEST(context.isNamed("Fluid"));
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(
        new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));

    boundingFromMapping->setMeshes(pSolidzMesh, pNastinMesh);
    boundingToMapping->setMeshes(pNastinMesh, pSolidzMesh);

    createNastinMesh3D(pNastinMesh);
    pNastinMesh->computeState();

    double safetyFactor = 0.0;

    ReceivedBoundingBox part(pSolidzMesh, safetyFactor);
    part.addM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicateBoundingBox();
    part.computeBoundingBox();
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI
