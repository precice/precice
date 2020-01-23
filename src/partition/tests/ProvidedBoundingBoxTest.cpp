#ifndef PRECICE_NO_MPI
#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicationFactory.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/SocketCommunication.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "m2n/PointToPointComFactory.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "partition/ProvidedBoundingBox.hpp"
#include "partition/ReceivedBoundingBox.hpp"
#include "partition/SharedPointer.hpp"
#include "testing/Fixtures.hpp"
#include "testing/Testing.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

using namespace precice;
using namespace partition;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ProvidedBoundingBoxTests)

void setupParallelEnvironment(m2n::PtrM2N m2n, const testing::TestContext &context)
{
  // Establish and configure the master-master connection
  if (context.isNamed("Fluid")) {
    BOOST_TEST(context.hasSize(3));
    if (context.isMaster()) {
      m2n->acceptMasterConnection("Fluid", "SolidMaster");
    }
  }

  if (context.isNamed("Solid")) {
    BOOST_TEST(context.hasSize(1));
    if (context.isMaster()) {
      m2n->requestMasterConnection("Fluid", "SolidMaster");
    }
  }

  if (context.isMaster()) {
      BOOST_TEST(m2n->isConnected());
  }
}

// create a communciator with two participants: each one has one master and one slave ranks
// only master-master and master-slave channels are created here
// Point to Point slave-slave channels must be created in the test.
void setupM2NBaseEnvironment(m2n::PtrM2N p2p, const testing::TestContext &context)
{
  // Both Solid and Fluid have 2 ranks
  BOOST_TEST(context.hasSize(2));

  // Establish and configure the master-master connection
  if (context.isNamed("Fluid")) {
    if (context.isMaster()) {
      p2p->acceptMasterConnection("FluidMaster", "SolidMaster");
    }
    utils::MasterSlave::configure(context.rank, context.size);
  }
  if (context.isNamed("Solid")) {
    if (context.isMaster()) { //Master Solid
      p2p->requestMasterConnection("FluidMaster", "SolidMaster");
    }
  }

  if (context.isMaster()) {
      BOOST_TEST(p2p->isConnected());
  }
}

void tearDownParallelEnvironment()
{
  mesh::Data::resetDataCount();
}

BOOST_AUTO_TEST_CASE(TestCommunicateBoundingBox2D)
{
  PRECICE_TEST("Fluid"_on(3_ranks).setupMasterSlaves(), "Solid"_on(1_rank), Require::Events);
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n, context);

  int  dimensions  = 2;
  bool flipNormals = true;

  if (context.isNamed("Fluid")) {

    mesh::PtrMesh pSolidMesh(new mesh::Mesh("SolidMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (context.isMaster()) { //Master
      Eigen::VectorXd position(dimensions);
      position << -1.0, 0.0;
      mesh::Vertex &v0 = pSolidMesh->createVertex(position);
      position << 1.0, 2.0;
      mesh::Vertex &v1 = pSolidMesh->createVertex(position);
      position << 5.0, 3.0;
      mesh::Vertex &v2 = pSolidMesh->createVertex(position);
      pSolidMesh->createEdge(v0, v1);
      pSolidMesh->createEdge(v1, v2);
    } else if (context.isRank(1)) { //Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.0, 3.5;
      mesh::Vertex &v3 = pSolidMesh->createVertex(position);
      position << 0.0, 4.5;
      mesh::Vertex &v4 = pSolidMesh->createVertex(position);
      pSolidMesh->createEdge(v3, v4);
    } else if (context.isRank(2)) { //Slave2
      Eigen::VectorXd position(dimensions);
      position << 2.5, 5.5;
      mesh::Vertex &v5 = pSolidMesh->createVertex(position);
      position << 4.5, 7.0;
      mesh::Vertex &v6 = pSolidMesh->createVertex(position);
      pSolidMesh->createEdge(v5, v6);
    }

    double safetyFactor = 0.0;
    bool   hasToSend    = true;
    pSolidMesh->computeState();

    ProvidedBoundingBox part(pSolidMesh, hasToSend, safetyFactor);
    part.addM2N(m2n);
    part.communicateBoundingBox();

  } else { //Solid
    BOOST_TEST(context.isNamed("Solid"));

    mesh::Mesh::BoundingBoxMap receivedGlobalBB;
    mesh::Mesh::BoundingBox    localBB;

    // we receive other participants communicator size
    int receivedFeedbackSize = 3;
    m2n->getMasterCommunication()->receive(receivedFeedbackSize, 0);

    for (int j = 0; j < dimensions; j++) {
      localBB.push_back(std::make_pair(-1, -1));
    }
    for (int i = 0; i < receivedFeedbackSize; i++) {
      receivedGlobalBB[i] = localBB;
    }

    // we receive golbal bounding box from othe participant!
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveBoundingBoxMap(receivedGlobalBB, 0);

    // check wether we have received the correct com size
    BOOST_TEST(receivedFeedbackSize == 3);

    //check the validity of received golbal bounding box (globalBB)
    BOOST_TEST(receivedGlobalBB[0][0].first == -1);
    BOOST_TEST(receivedGlobalBB[0][0].second == 5);
    BOOST_TEST(receivedGlobalBB[0][1].first == 0);
    BOOST_TEST(receivedGlobalBB[0][1].second == 3);
    BOOST_TEST(receivedGlobalBB[1][0].first == 0);
    BOOST_TEST(receivedGlobalBB[1][0].second == 1);
    BOOST_TEST(receivedGlobalBB[1][1].first == 3.5);
    BOOST_TEST(receivedGlobalBB[1][1].second == 4.5);
    BOOST_TEST(receivedGlobalBB[2][0].first == 2.5);
    BOOST_TEST(receivedGlobalBB[2][0].second == 4.5);
    BOOST_TEST(receivedGlobalBB[2][1].first == 5.5);
    BOOST_TEST(receivedGlobalBB[2][1].second == 7.0);
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestCommunicateBoundingBox3D)
{
  PRECICE_TEST("Fluid"_on(3_ranks).setupMasterSlaves(), "Solid"_on(1_rank), Require::Events);
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n, context);

  int  dimensions  = 3;
  bool flipNormals = true;

  if (context.isNamed("Fluid")) { //NASTIN

    mesh::PtrMesh pSolidMesh(new mesh::Mesh("SolidMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (context.isMaster()) { //Master
      Eigen::VectorXd position(dimensions);
      position << -1.0, 0.0, -1.0;
      mesh::Vertex &v0 = pSolidMesh->createVertex(position);
      position << 1.0, 2.0, 1.0;
      mesh::Vertex &v1 = pSolidMesh->createVertex(position);
      position << 5.0, 3.0, 5.0;
      mesh::Vertex &v2 = pSolidMesh->createVertex(position);
      pSolidMesh->createEdge(v0, v1);
      pSolidMesh->createEdge(v1, v2);
    }

    else if (context.isRank(1)) { //Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.0, 3.5, 1.0;
      mesh::Vertex &v3 = pSolidMesh->createVertex(position);
      position << 0.0, 4.5, 0.0;
      mesh::Vertex &v4 = pSolidMesh->createVertex(position);
      pSolidMesh->createEdge(v3, v4);
    } else if (context.isRank(2)) { //Slave2
      Eigen::VectorXd position(dimensions);
      position << 2.5, 5.5, 2.5;
      mesh::Vertex &v5 = pSolidMesh->createVertex(position);
      position << 4.5, 7.0, 4.5;
      mesh::Vertex &v6 = pSolidMesh->createVertex(position);
      pSolidMesh->createEdge(v5, v6);
    }

    double safetyFactor = 0.0;
    bool   hasToSend    = true;
    pSolidMesh->computeState();

    ProvidedBoundingBox part(pSolidMesh, hasToSend, safetyFactor);
    part.addM2N(m2n);
    part.communicateBoundingBox();

  } else { //Solid

    mesh::Mesh::BoundingBoxMap receivedGlobalBB;
    mesh::Mesh::BoundingBox    localBB;

    // we receive other participants communicator size
    int remoteParComSize = 3;
    m2n->getMasterCommunication()->receive(remoteParComSize, 0);

    for (int j = 0; j < dimensions; j++) {
      localBB.push_back(std::make_pair(-1, -1));
    }
    for (int i = 0; i < remoteParComSize; i++) {
      receivedGlobalBB[i] = localBB;
    }

    // we receive golbal bounding box from othe participant!
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveBoundingBoxMap(receivedGlobalBB, 0);

    // check wether we have received the correct com size
    BOOST_TEST(remoteParComSize == 3);

    //check the validity of received golbal bounding box (globalBB)
    BOOST_TEST(receivedGlobalBB[0][0].first == -1);
    BOOST_TEST(receivedGlobalBB[0][0].second == 5);
    BOOST_TEST(receivedGlobalBB[0][1].first == 0);
    BOOST_TEST(receivedGlobalBB[0][1].second == 3);
    BOOST_TEST(receivedGlobalBB[0][2].first == -1);
    BOOST_TEST(receivedGlobalBB[0][2].second == 5);
    BOOST_TEST(receivedGlobalBB[1][0].first == 0);
    BOOST_TEST(receivedGlobalBB[1][0].second == 1);
    BOOST_TEST(receivedGlobalBB[1][1].first == 3.5);
    BOOST_TEST(receivedGlobalBB[1][1].second == 4.5);
    BOOST_TEST(receivedGlobalBB[1][2].first == 0);
    BOOST_TEST(receivedGlobalBB[1][2].second == 1);
    BOOST_TEST(receivedGlobalBB[2][0].first == 2.5);
    BOOST_TEST(receivedGlobalBB[2][0].second == 4.5);
    BOOST_TEST(receivedGlobalBB[2][1].first == 5.5);
    BOOST_TEST(receivedGlobalBB[2][1].second == 7.0);
    BOOST_TEST(receivedGlobalBB[2][2].first == 2.5);
    BOOST_TEST(receivedGlobalBB[2][2].second == 4.5);
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestComputeBoundingBox)
{
  PRECICE_TEST("Fluid"_on(3_ranks).setupMasterSlaves(), "Solid"_on(1_rank), Require::Events);
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n, context);

  int  dimensions  = 2;
  bool flipNormals = true;

  if (context.isNamed("Solid")) {
    std::vector<int> connectedRanks = {0, 1, 2};
    m2n->getMasterCommunication()->send(connectedRanks, 0);

    // construct connection map
    std::map<int, std::vector<int>> sendConnectionMap;
    sendConnectionMap[0].push_back(1);
    sendConnectionMap[0].push_back(2);
    sendConnectionMap[1].push_back(0);
    sendConnectionMap[1].push_back(2);
    sendConnectionMap[2].push_back(0);
    sendConnectionMap[2].push_back(1);

    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendConnectionMap(sendConnectionMap, 0);
  } else { //Fluid

    mesh::PtrMesh pSolidMesh(new mesh::Mesh("SolidMesh", dimensions, flipNormals, testing::nextMeshID()));

    double safetyFactor = 0.0;
    bool   hasToSend    = true;
    pSolidMesh->computeState();

    ProvidedBoundingBox part(pSolidMesh, hasToSend, safetyFactor);
    part.addM2N(m2n);
    part.computeBoundingBox();

    if (context.isMaster()) { //Master
      BOOST_TEST(pSolidMesh->getConnectedRanks().size() == 2);
      BOOST_TEST(pSolidMesh->getConnectedRanks()[0] == 1);
      BOOST_TEST(pSolidMesh->getConnectedRanks()[1] == 2);
    } else if (context.isRank(1)) { //Slave1
      BOOST_TEST(pSolidMesh->getConnectedRanks().size() == 2);
      BOOST_TEST(pSolidMesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(pSolidMesh->getConnectedRanks()[1] == 2);
    } else if (context.isRank(2)) { //Slave2
      BOOST_TEST(pSolidMesh->getConnectedRanks().size() == 2);
      BOOST_TEST(pSolidMesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(pSolidMesh->getConnectedRanks()[1] == 1);
    }
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestCommunicateLocalMeshPartitions)
{
  PRECICE_TEST("Solid"_on(2_ranks).setupMasterSlaves(), "Fluid"_on(2_ranks).setupMasterSlaves(), Require::Events);
  //mesh creation
  int           dimensions   = 2;
  bool          flipNormals  = true;
  double        safetyFactor = 0.1;
  bool          hasToSend    = true;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));

  // create second communicator for m2n mesh and communciation map exchange
  com::PtrCommunication                     participantsCom       = com::PtrCommunication(new com::SocketCommunication());
  com::PtrCommunicationFactory              participantComFactory = com::PtrCommunicationFactory(new com::SocketCommunicationFactory);
  m2n::DistributedComFactory::SharedPointer distributionFactory   = m2n::DistributedComFactory::SharedPointer(new m2n::PointToPointComFactory(participantComFactory));
  m2n::PtrM2N                               p2p                   = m2n::PtrM2N(new m2n::M2N(participantsCom, distributionFactory));

  setupM2NBaseEnvironment(p2p, context);

  if (context.isNamed("Solid")) {
    if (context.isMaster()) {
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

    } else {
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
    }
  } else {
    // Fluid
    if (context.isMaster()) {
      mesh->getConnectedRanks().push_back(0);
    } else if (context.isRank(1)) {
      mesh->getConnectedRanks().push_back(1);
    }
  }

  mesh->computeState();

  if (context.isNamed("Solid")) {
    p2p->createDistributedCommunication(mesh);
    ProvidedBoundingBox part(mesh, hasToSend, safetyFactor);
    p2p->acceptSlavesPreConnection("SolidSlaves", "FluidSlaves");
    part.addM2N(p2p);

    part.communicate();
  } else { // Fluid
    p2p->createDistributedCommunication(mesh);
    ReceivedBoundingBox part(mesh, safetyFactor);
    p2p->requestSlavesPreConnection("SolidSlaves", "FluidSlaves");
    part.addM2N(p2p);

    part.communicate();

    BOOST_TEST(mesh->vertices().size() == 4);

    if (context.isMaster()) {
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

BOOST_AUTO_TEST_CASE(TestCompute2D)
{
  PRECICE_TEST("Solid"_on(2_ranks).setupMasterSlaves(), "Fluid"_on(2_ranks).setupMasterSlaves(), Require::Events);
  //mesh creation
  int           dimensions   = 2;
  bool          flipNormals  = true;
  double        safetyFactor = 0;
  bool          hasToSend    = true;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh receivedMesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));

  if (context.isNamed("Solid")) {
    if (context.isMaster()) {
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
    } else {
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
    }
  } else {
    // Fluid
    if (context.isMaster()) {
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
    } else {
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
    }
  }

  mesh->computeState();

  // create the communicator for m2n mesh and communciation map exchange
  com::PtrCommunication                     participantsCom       = com::PtrCommunication(new com::SocketCommunication());
  com::PtrCommunicationFactory              participantComFactory = com::PtrCommunicationFactory(new com::SocketCommunicationFactory);
  m2n::DistributedComFactory::SharedPointer distributionFactory   = m2n::DistributedComFactory::SharedPointer(new m2n::PointToPointComFactory(participantComFactory));
  m2n::PtrM2N                               p2p                   = m2n::PtrM2N(new m2n::M2N(participantsCom, distributionFactory));

  setupM2NBaseEnvironment(p2p, context);

  if (context.isNamed("Solid")) {
    p2p->createDistributedCommunication(mesh);
    ProvidedBoundingBox part(mesh, hasToSend, safetyFactor);
    part.addM2N(p2p);

    part.communicateBoundingBox();
    part.computeBoundingBox();

    if (context.isMaster()) {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(mesh->getConnectedRanks()[1] == 1);
    } else {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(mesh->getConnectedRanks()[1] == 1);
    }

    p2p->acceptSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();

    if (context.isMaster()) {
      BOOST_TEST(mesh->getCommunicationMap()[0][0] == 0);
      BOOST_TEST(mesh->getCommunicationMap()[0][1] == 1);
    } else {
      BOOST_TEST(mesh->getCommunicationMap()[0][0] == 0);
      BOOST_TEST(mesh->getCommunicationMap()[1][0] == 1);
      BOOST_TEST(mesh->getCommunicationMap()[1][1] == 2);
    }

  } else {
    // Fluid
    p2p->createDistributedCommunication(receivedMesh);
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping   = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(receivedMesh, mesh);
    boundingToMapping->setMeshes(mesh, receivedMesh);

    ReceivedBoundingBox part(receivedMesh, safetyFactor);

    part.addM2N(p2p);

    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);

    part.communicateBoundingBox();
    part.computeBoundingBox();

    p2p->requestSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestCompute3D)
{
  PRECICE_TEST("Solid"_on(2_ranks).setupMasterSlaves(), "Fluid"_on(2_ranks).setupMasterSlaves(), Require::Events);
  //mesh creation
  int           dimensions   = 3;
  bool          flipNormals  = true;
  double        safetyFactor = 0.0;
  bool          hasToSend    = true;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh receivedMesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));

  if (context.isNamed("Solid")) {
    if (context.isMaster()) {
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
    } else {
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
    }
  } else {
    // Fluid
    if (context.isMaster()) {
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
    } else {
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
    }
  }

  mesh->computeState();

  // create the communicator for m2n mesh and communciation map exchange
  com::PtrCommunication                     participantsCom       = com::PtrCommunication(new com::SocketCommunication());
  com::PtrCommunicationFactory              participantComFactory = com::PtrCommunicationFactory(new com::SocketCommunicationFactory);
  m2n::DistributedComFactory::SharedPointer distributionFactory   = m2n::DistributedComFactory::SharedPointer(new m2n::PointToPointComFactory(participantComFactory));
  m2n::PtrM2N                               p2p                   = m2n::PtrM2N(new m2n::M2N(participantsCom, distributionFactory));

  setupM2NBaseEnvironment(p2p, context);

  if (context.isNamed("Solid")) {
    p2p->createDistributedCommunication(mesh);
    ProvidedBoundingBox part(mesh, hasToSend, safetyFactor);
    part.addM2N(p2p);

    part.communicateBoundingBox();
    part.computeBoundingBox();

    if (context.isMaster()) {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(mesh->getConnectedRanks()[1] == 1);
    } else {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks()[0] == 0);
      BOOST_TEST(mesh->getConnectedRanks()[1] == 1);
    }

    p2p->acceptSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();

    if (context.isMaster()) {
      BOOST_TEST(mesh->getCommunicationMap()[0][0] == 0);
      BOOST_TEST(mesh->getCommunicationMap()[0][1] == 1);
    } else {
      BOOST_TEST(mesh->getCommunicationMap()[0][0] == 0);
      BOOST_TEST(mesh->getCommunicationMap()[1][0] == 1);
      BOOST_TEST(mesh->getCommunicationMap()[1][1] == 2);
    }
  } else {
    // Fluid
    p2p->createDistributedCommunication(receivedMesh);
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping   = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(receivedMesh, mesh);
    boundingToMapping->setMeshes(mesh, receivedMesh);

    ReceivedBoundingBox part(receivedMesh, safetyFactor);
    part.addM2N(p2p);

    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicateBoundingBox();
    part.computeBoundingBox();

    p2p->requestSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI
