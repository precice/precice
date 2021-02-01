#ifndef PRECICE_NO_MPI
#include <Eigen/Core>
#include <algorithm>
#include <deque>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "com/CommunicateBoundingBox.hpp"
#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "m2n/M2N.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/SharedPointer.hpp"
#include "math/constants.hpp"
#include "mesh/BoundingBox.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "partition/Partition.hpp"
#include "partition/ProvidedPartition.hpp"
#include "partition/ReceivedPartition.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mesh {
class Edge;
class Vertex;
} // namespace mesh
} // namespace precice

using namespace precice;
using namespace partition;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ProvidedPartitionTests)

void tearDownParallelEnvironment()
{
  mesh::Data::resetDataCount();
}

BOOST_AUTO_TEST_CASE(TestGatherAndCommunicate2D)
{
  PRECICE_TEST("NASTIN"_on(1_rank), "SOLIDZ"_on(3_ranks).setupMasterSlaves(), Require::Events);
  auto m2n = context.connectMasters("NASTIN", "SOLIDZ");

  int  dimensions  = 2;
  bool flipNormals = false;

  if (context.isNamed("NASTIN")) { //NASTIN
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    double safetyFactor = 0.1;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_SLAVES, safetyFactor);
    part.addM2N(m2n);
    part.communicate();

    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 4);

    for (int i = 0; i < 6; i++) {
      BOOST_TEST(pSolidzMesh->vertices().at(i).getGlobalIndex() == i);
    }
  } else { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (context.isMaster()) { //Master
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0;
      mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
      position << 0.0, 1.5;
      mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v1, v2);
    } else if (context.isRank(1)) { //Slave1
    } else if (context.isRank(2)) { //Slave2
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
    pSolidzMesh->computeBoundingBox();
    pSolidzMesh->computeState();

    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    BOOST_REQUIRE(pSolidzMesh->getVertexOffsets().size() == 3);
    if (context.isMaster()) { //master
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(1) == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(2) == 6);
    } else if (context.isRank(1)) { //Slave1
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(1) == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(2) == 6);
    } else {
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(0) == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(1) == 2);
      BOOST_TEST(pSolidzMesh->getVertexOffsets().at(2) == 6);
    }
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestGatherAndCommunicate3D)
{
  PRECICE_TEST("NASTIN"_on(1_rank), "SOLIDZ"_on(3_ranks).setupMasterSlaves(), Require::Events);
  auto m2n = context.connectMasters("NASTIN", "SOLIDZ");

  int  dimensions  = 3;
  bool flipNormals = false;

  if (context.isNamed("NASTIN")) { //NASTIN
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    double safetyFactor = 0.1;

    ReceivedPartition part(pSolidzMesh, ReceivedPartition::ON_SLAVES, safetyFactor);
    part.addM2N(m2n);
    part.communicate();

    BOOST_TEST(pSolidzMesh->vertices().size() == 6);
    BOOST_TEST(pSolidzMesh->edges().size() == 6);
    BOOST_TEST(pSolidzMesh->triangles().size() == 2);

    for (int i = 0; i < 6; i++) {
      BOOST_TEST(pSolidzMesh->vertices().at(i).getGlobalIndex() == i);
    }
  } else { //SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (context.isMaster()) { //Master
      Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0, 0.0;
      mesh::Vertex &v1 = pSolidzMesh->createVertex(position);
      position << 0.0, 1.5, 1.0;
      mesh::Vertex &v2 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v1, v2);
    } else if (context.isRank(1)) { //Slave1
    } else if (context.isRank(2)) { //Slave2
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
    pSolidzMesh->computeBoundingBox();
    pSolidzMesh->computeState();

    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.communicate();
    part.compute();

    BOOST_TEST(pSolidzMesh->getGlobalNumberOfVertices() == 6);
    const auto &vertices           = pSolidzMesh->vertices();
    const auto &vertexOffsets      = pSolidzMesh->getVertexOffsets();
    const auto &vertexDistribution = pSolidzMesh->getVertexDistribution();

    if (context.isMaster()) { //master
      BOOST_REQUIRE(vertexOffsets.size() == 3);
      BOOST_TEST(vertexOffsets.at(0) == 2);
      BOOST_TEST(vertexOffsets.at(1) == 2);
      BOOST_TEST(vertexOffsets.at(2) == 6);

      BOOST_REQUIRE(vertices.size() == 2);
      BOOST_TEST(vertices.at(0).getGlobalIndex() == 0);
      BOOST_TEST(vertices.at(1).getGlobalIndex() == 1);
      BOOST_TEST(vertices.at(0).isOwner() == true);
      BOOST_TEST(vertices.at(1).isOwner() == true);

      BOOST_REQUIRE((vertexDistribution.size()) == 3);
      BOOST_TEST(vertexDistribution.at(0).size() == 2);
      BOOST_TEST(vertexDistribution.at(1).size() == 0);
      BOOST_TEST(vertexDistribution.at(2).size() == 4);
      BOOST_TEST(vertexDistribution.at(0).at(0) == 0);
      BOOST_TEST(vertexDistribution.at(0).at(1) == 1);
      BOOST_TEST(vertexDistribution.at(2).at(0) == 2);
      BOOST_TEST(vertexDistribution.at(2).at(1) == 3);
      BOOST_TEST(vertexDistribution.at(2).at(2) == 4);
      BOOST_TEST(vertexDistribution.at(2).at(3) == 5);
    } else if (context.isRank(1)) { //Slave1
      BOOST_REQUIRE(vertexOffsets.size() == 3);
      BOOST_TEST(vertexOffsets.at(0) == 2);
      BOOST_TEST(vertexOffsets.at(1) == 2);
      BOOST_TEST(vertexOffsets.at(2) == 6);
    } else if (context.isRank(2)) { //Slave2
      BOOST_REQUIRE(vertexOffsets.size() == 3);
      BOOST_TEST(vertexOffsets.at(0) == 2);
      BOOST_TEST(vertexOffsets.at(1) == 2);
      BOOST_TEST(vertexOffsets.at(2) == 6);

      BOOST_REQUIRE(vertices.size() == 4);
      BOOST_TEST(vertices.at(0).getGlobalIndex() == 2);
      BOOST_TEST(vertices.at(1).getGlobalIndex() == 3);
      BOOST_TEST(vertices.at(2).getGlobalIndex() == 4);
      BOOST_TEST(vertices.at(3).getGlobalIndex() == 5);
      BOOST_TEST(vertices.at(0).isOwner() == true);
      BOOST_TEST(vertices.at(1).isOwner() == true);
      BOOST_TEST(vertices.at(2).isOwner() == true);
      BOOST_TEST(vertices.at(3).isOwner() == true);
    }
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestOnlyDistribution2D)
{
  PRECICE_TEST("NASTIN"_on(4_ranks).setupMasterSlaves(), Require::Events);
  // Create mesh object
  std::string   meshName("MyMesh");
  int           dim         = 2;
  bool          flipNormals = false; // The normals of triangles, edges, vertices
  mesh::PtrMesh pMesh(new mesh::Mesh(meshName, dim, flipNormals, testing::nextMeshID()));

  if (context.isMaster()) { //Master
    Eigen::VectorXd position(dim);
    position << 0.0, 0.0;
    pMesh->createVertex(position);
    position << 1.0, 0.0;
    pMesh->createVertex(position);
  } else if (context.isRank(1)) { //Slave1
    Eigen::VectorXd position(dim);
    position << 2.0, 0.0;
    pMesh->createVertex(position);
  } else if (context.isRank(2)) { //Slave2
  } else if (context.isRank(3)) { //Slave3
    Eigen::VectorXd position(dim);
    position << 3.0, 0.0;
    pMesh->createVertex(position);
    position << 4.0, 0.0;
    pMesh->createVertex(position);
  }
  pMesh->computeState();
  pMesh->computeBoundingBox();

  ProvidedPartition part(pMesh);
  part.communicate();
  part.compute();

  BOOST_TEST_CONTEXT(*pMesh)
  {
    if (context.isMaster()) { //Master
      BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 5);
      BOOST_TEST_REQUIRE(pMesh->getVertexOffsets().size() == 4);
      BOOST_TEST(pMesh->getVertexOffsets().at(0) == 2);
      BOOST_TEST(pMesh->getVertexOffsets().at(1) == 3);
      BOOST_TEST(pMesh->getVertexOffsets().at(2) == 3);
      BOOST_TEST(pMesh->getVertexOffsets().at(3) == 5);
      BOOST_TEST(pMesh->vertices().at(0).getGlobalIndex() == 0);
      BOOST_TEST(pMesh->vertices().at(1).getGlobalIndex() == 1);
      BOOST_TEST(pMesh->vertices().at(0).isOwner() == true);
      BOOST_TEST(pMesh->vertices().at(1).isOwner() == true);
    } else if (context.isRank(1)) { //Slave1
      BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 5);
      BOOST_TEST_REQUIRE(pMesh->getVertexOffsets().size() == 4);
      BOOST_TEST(pMesh->getVertexOffsets().at(0) == 2);
      BOOST_TEST(pMesh->getVertexOffsets().at(1) == 3);
      BOOST_TEST(pMesh->getVertexOffsets().at(2) == 3);
      BOOST_TEST(pMesh->getVertexOffsets().at(3) == 5);
      BOOST_TEST(pMesh->vertices().at(0).getGlobalIndex() == 2);
      BOOST_TEST(pMesh->vertices().at(0).isOwner() == true);
    } else if (context.isRank(2)) { //Slave2
      BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 5);
      BOOST_TEST_REQUIRE(pMesh->getVertexOffsets().size() == 4);
      BOOST_TEST(pMesh->getVertexOffsets().at(0) == 2);
      BOOST_TEST(pMesh->getVertexOffsets().at(1) == 3);
      BOOST_TEST(pMesh->getVertexOffsets().at(2) == 3);
      BOOST_TEST(pMesh->getVertexOffsets().at(3) == 5);
    } else if (context.isRank(3)) { //Slave3
      BOOST_TEST(pMesh->getGlobalNumberOfVertices() == 5);
      BOOST_TEST_REQUIRE(pMesh->getVertexOffsets().size() == 4);
      BOOST_TEST(pMesh->getVertexOffsets().at(0) == 2);
      BOOST_TEST(pMesh->getVertexOffsets().at(1) == 3);
      BOOST_TEST(pMesh->getVertexOffsets().at(2) == 3);
      BOOST_TEST(pMesh->getVertexOffsets().at(3) == 5);
      BOOST_TEST(pMesh->vertices().at(0).getGlobalIndex() == 3);
      BOOST_TEST(pMesh->vertices().at(1).getGlobalIndex() == 4);
      BOOST_TEST(pMesh->vertices().at(0).isOwner() == true);
      BOOST_TEST(pMesh->vertices().at(1).isOwner() == true);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestCompareBoundingBoxes2D)
{
  PRECICE_TEST("SOLIDZ"_on(3_ranks).setupMasterSlaves(), "NASTIN"_on(1_rank), Require::Events);
  testing::ConnectionOptions options;
  options.useOnlyMasterCom = false;
  options.useTwoLevelInit  = true;
  auto m2n                 = context.connectMasters("NASTIN", "SOLIDZ", options);

  int  dimensions  = 2;
  bool flipNormals = true;

  if (context.isNamed("SOLIDZ")) { //SOLIDZ

    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (context.isMaster()) { //Master
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

    else if (context.isRank(1)) { //Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.0, 3.5;
      mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5;
      mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3, v4);
    } else if (context.isRank(2)) { //Slave2
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

    if (context.isMaster()) { //Master
      BOOST_TEST(pSolidzMesh->getConnectedRanks().size() == 2);
      BOOST_TEST(pSolidzMesh->getConnectedRanks().at(0) == 1);
      BOOST_TEST(pSolidzMesh->getConnectedRanks().at(1) == 2);
    } else if (context.isRank(1)) { //Slave1
      BOOST_TEST(pSolidzMesh->getConnectedRanks().size() == 2);
      BOOST_TEST(pSolidzMesh->getConnectedRanks().at(0) == 0);
      BOOST_TEST(pSolidzMesh->getConnectedRanks().at(1) == 2);
    } else if (context.isRank(2)) { //Slave2
      BOOST_TEST(pSolidzMesh->getConnectedRanks().size() == 2);
      BOOST_TEST(pSolidzMesh->getConnectedRanks().at(0) == 0);
      BOOST_TEST(pSolidzMesh->getConnectedRanks().at(1) == 1);
    }

  } else { //NASTIN
    BOOST_TEST(context.isNamed("NASTIN"));

    mesh::Mesh::BoundingBoxMap receivedGlobalBB;
    mesh::BoundingBox          localBB{dimensions};

    mesh::Mesh::BoundingBoxMap compareBB;
    compareBB.emplace(0, mesh::BoundingBox({-1, 5, 0, 3}));
    compareBB.emplace(1, mesh::BoundingBox({0, 1, 3.5, 4.5}));
    compareBB.emplace(2, mesh::BoundingBox({2.5, 4.5, 5.5, 7.0}));

    // we receive other participants communicator size
    int receivedFeedbackSize = 3;
    m2n->getMasterCommunication()->receive(receivedFeedbackSize, 0);

    for (int i = 0; i < receivedFeedbackSize; i++) {
      receivedGlobalBB.emplace(i, localBB);
    }

    // we receive golbal bounding box from othe participant!
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveBoundingBoxMap(receivedGlobalBB, 0);
    // check wether we have received the correct com size
    BOOST_TEST(receivedFeedbackSize == 3);

    //check the validity of received golbal bounding box (globalBB)
    BOOST_TEST(receivedGlobalBB.at(0) == compareBB.at(0));
    BOOST_TEST(receivedGlobalBB.at(1) == compareBB.at(1));
    BOOST_TEST(receivedGlobalBB.at(2) == compareBB.at(2));

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

BOOST_AUTO_TEST_CASE(TestSendBoundingBoxes3D)
{
  PRECICE_TEST("SOLIDZ"_on(3_ranks).setupMasterSlaves(), "NASTIN"_on(1_rank), Require::Events);
  testing::ConnectionOptions options;
  options.useOnlyMasterCom = false;
  options.useTwoLevelInit  = true;
  auto m2n                 = context.connectMasters("NASTIN", "SOLIDZ", options);

  int  dimensions  = 3;
  bool flipNormals = true;

  if (context.isNamed("SOLIDZ")) { //SOLIDZ

    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals, testing::nextMeshID()));

    if (context.isMaster()) { //Master
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

    else if (context.isRank(1)) { //Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.0, 3.5, 1.0;
      mesh::Vertex &v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5, 0.0;
      mesh::Vertex &v4 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3, v4);
    } else if (context.isRank(2)) { //Slave2
      Eigen::VectorXd position(dimensions);
      position << 2.5, 5.5, 2.5;
      mesh::Vertex &v5 = pSolidzMesh->createVertex(position);
      position << 4.5, 7.0, 4.5;
      mesh::Vertex &v6 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v5, v6);
    }
    pSolidzMesh->computeState();
    pSolidzMesh->computeBoundingBox();

    ProvidedPartition part(pSolidzMesh);
    part.addM2N(m2n);
    part.compareBoundingBoxes();

  } else { //NASTIN
    BOOST_TEST(context.isNamed("NASTIN"));

    mesh::Mesh::BoundingBoxMap receivedGlobalBB;
    mesh::BoundingBox          localBB{dimensions};

    mesh::Mesh::BoundingBoxMap compareBB;
    compareBB.emplace(0, mesh::BoundingBox({-1, 5, 0, 3, -1, 5}));
    compareBB.emplace(1, mesh::BoundingBox({0, 1, 3.5, 4.5, 0, 1}));
    compareBB.emplace(2, mesh::BoundingBox({2.5, 4.5, 5.5, 7.0, 2.5, 4.5}));

    // we receive other participants communicator size
    int remoteParComSize = 3;
    m2n->getMasterCommunication()->receive(remoteParComSize, 0);

    for (int i = 0; i < remoteParComSize; i++) {
      receivedGlobalBB.emplace(i, localBB);
    }

    // we receive golbal bounding box from othe participant!
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveBoundingBoxMap(receivedGlobalBB, 0);

    // check wether we have received the correct com size
    BOOST_TEST(remoteParComSize == 3);

    //check the validity of received golbal bounding box (globalBB)
    BOOST_TEST(receivedGlobalBB.at(0) == compareBB.at(0));
    BOOST_TEST(receivedGlobalBB.at(1) == compareBB.at(1));
    BOOST_TEST(receivedGlobalBB.at(2) == compareBB.at(2));

    //send empty dummy list of connected ranks as feedback
    std::vector<int> connectedRanksList;
    m2n->getMasterCommunication()->send(connectedRanksList, 0);
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
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));

  testing::ConnectionOptions options;
  options.useOnlyMasterCom = false;
  options.useTwoLevelInit  = true;
  options.type             = testing::ConnectionType::PointToPoint;
  auto m2n                 = context.connectMasters("Fluid", "Solid", options);

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
    BOOST_TEST(context.isNamed("Fluid"));
    if (context.isMaster()) {
      mesh->getConnectedRanks().push_back(0);
    } else {
      mesh->getConnectedRanks().push_back(1);
    }
  }
  mesh->computeBoundingBox();
  mesh->computeState();

  if (context.isNamed("Solid")) {
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

    if (context.isMaster()) {
      BOOST_TEST(mesh->vertices().at(0).getCoords()(0) == 0.5);
      BOOST_TEST(mesh->vertices().at(0).getCoords()(1) == 0.0);
      BOOST_TEST(mesh->vertices().at(1).getCoords()(0) == 1.5);
      BOOST_TEST(mesh->vertices().at(1).getCoords()(1) == 0.0);
      BOOST_TEST(mesh->vertices().at(2).getCoords()(0) == 2.0);
      BOOST_TEST(mesh->vertices().at(2).getCoords()(1) == 1.0);
      BOOST_TEST(mesh->vertices().at(3).getCoords()(0) == 0.5);
      BOOST_TEST(mesh->vertices().at(3).getCoords()(1) == 1.0);
    } else {
      BOOST_TEST(mesh->vertices().at(0).getCoords()(0) == 2.5);
      BOOST_TEST(mesh->vertices().at(0).getCoords()(1) == 0.0);
      BOOST_TEST(mesh->vertices().at(1).getCoords()(0) == 3.5);
      BOOST_TEST(mesh->vertices().at(1).getCoords()(1) == 0.0);
      BOOST_TEST(mesh->vertices().at(2).getCoords()(0) == 3.5);
      BOOST_TEST(mesh->vertices().at(2).getCoords()(1) == 1.0);
      BOOST_TEST(mesh->vertices().at(3).getCoords()(0) == 2.0);
      BOOST_TEST(mesh->vertices().at(3).getCoords()(1) == 1.0);
    }
  }

  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestTwoLevelRepartitioning2D)
{
  PRECICE_TEST("Solid"_on(2_ranks).setupMasterSlaves(), "Fluid"_on(2_ranks).setupMasterSlaves(), Require::Events);
  //mesh creation
  int           dimensions   = 2;
  bool          flipNormals  = true;
  double        safetyFactor = 0;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh receivedMesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));

  testing::ConnectionOptions options;
  options.useOnlyMasterCom = false;
  options.useTwoLevelInit  = true;
  options.type             = testing::ConnectionType::PointToPoint;
  auto m2n                 = context.connectMasters("Fluid", "Solid", options);

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
    BOOST_TEST(context.isNamed("Fluid"));
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
  mesh->computeBoundingBox();

  if (context.isNamed("Solid")) {
    m2n->createDistributedCommunication(mesh);
    ProvidedPartition part(mesh);
    part.addM2N(m2n);

    part.compareBoundingBoxes();

    if (context.isMaster()) {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks().at(0) == 0);
      BOOST_TEST(mesh->getConnectedRanks().at(1) == 1);
    } else {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks().at(0) == 0);
      BOOST_TEST(mesh->getConnectedRanks().at(1) == 1);
    }

    m2n->acceptSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();

    if (context.isMaster()) {
      BOOST_TEST(mesh->getCommunicationMap().at(0).at(0) == 0);
      BOOST_TEST(mesh->getCommunicationMap().at(0).at(1) == 1);
    } else {
      BOOST_TEST(mesh->getCommunicationMap().at(0).at(0) == 0);
      BOOST_TEST(mesh->getCommunicationMap().at(1).at(0) == 1);
      BOOST_TEST(mesh->getCommunicationMap().at(1).at(1) == 2);
    }
  } else {
    m2n->createDistributedCommunication(receivedMesh);
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping   = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(receivedMesh, mesh);
    boundingToMapping->setMeshes(mesh, receivedMesh);

    ReceivedPartition part(receivedMesh, ReceivedPartition::ON_SLAVES, safetyFactor);

    part.addM2N(m2n);

    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);

    part.compareBoundingBoxes();

    m2n->requestSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();
  }
  tearDownParallelEnvironment();
}

BOOST_AUTO_TEST_CASE(TestTwoLevelRepartitioning3D)
{
  PRECICE_TEST("Solid"_on(2_ranks).setupMasterSlaves(), "Fluid"_on(2_ranks).setupMasterSlaves(), Require::Events);

  //mesh creation
  int           dimensions   = 3;
  bool          flipNormals  = true;
  double        safetyFactor = 0.0;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));
  mesh::PtrMesh receivedMesh(new mesh::Mesh("mesh", dimensions, flipNormals, testing::nextMeshID()));

  // create the communicator for m2n mesh and communciation map exchange
  testing::ConnectionOptions options;
  options.useOnlyMasterCom = false;
  options.useTwoLevelInit  = true;
  options.type             = testing::ConnectionType::PointToPoint;
  auto m2n                 = context.connectMasters("Fluid", "Solid", options);

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
  mesh->computeBoundingBox();

  if (context.isNamed("Solid")) {
    m2n->createDistributedCommunication(mesh);
    ProvidedPartition part(mesh);
    part.addM2N(m2n);

    part.compareBoundingBoxes();

    if (context.isMaster()) {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks().at(0) == 0);
      BOOST_TEST(mesh->getConnectedRanks().at(1) == 1);
    } else {
      BOOST_TEST(mesh->getConnectedRanks().size() == 2);
      BOOST_TEST(mesh->getConnectedRanks().at(0) == 0);
      BOOST_TEST(mesh->getConnectedRanks().at(1) == 1);
    }

    m2n->acceptSlavesPreConnection("FluidSlaves", "SolidSlaves");

    part.communicate();
    part.compute();

    if (context.isMaster()) {
      BOOST_TEST(mesh->getCommunicationMap().at(0).at(0) == 0);
      BOOST_TEST(mesh->getCommunicationMap().at(0).at(1) == 1);
    } else {
      BOOST_TEST(mesh->getCommunicationMap().at(0).at(0) == 0);
      BOOST_TEST(mesh->getCommunicationMap().at(1).at(0) == 1);
      BOOST_TEST(mesh->getCommunicationMap().at(1).at(1) == 2);
    }
  } else {
    m2n->createDistributedCommunication(receivedMesh);
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping   = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(receivedMesh, mesh);
    boundingToMapping->setMeshes(mesh, receivedMesh);

    ReceivedPartition part(receivedMesh, ReceivedPartition::ON_SLAVES, safetyFactor);
    part.addM2N(m2n);

    part.addFromMapping(boundingFromMapping);
    part.addToMapping(boundingToMapping);

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
