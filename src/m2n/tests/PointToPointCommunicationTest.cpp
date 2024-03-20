#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <algorithm>
#include <map>
#include <memory>
#include <vector>
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/SharedPointer.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "m2n/PointToPointCommunication.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/IntraComm.hpp"

namespace precice {
namespace mesh {
class Vertex;
} // namespace mesh
namespace utils {
class Parallel;
} // namespace utils
} // namespace precice

using precice::testing::TestContext;
using precice::utils::IntraComm;
using precice::utils::Parallel;

using std::vector;

using namespace precice;
using namespace m2n;

BOOST_AUTO_TEST_SUITE(M2NTests)

void process(vector<double> &data)
{
  for (auto &elem : data) {
    elem += IntraComm::getRank() + 1;
  }
}

void runP2PComTest1(const TestContext &context, com::PtrCommunicationFactory cf)
{
  BOOST_TEST(context.hasSize(2));

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));

  m2n::PointToPointCommunication c(cf, mesh);

  vector<double> data;
  vector<double> expectedData;

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      mesh->setGlobalNumberOfVertices(10);

      mesh->setVertexDistribution({{0, {0, 1, 3, 5, 7}}, {1, {1, 2, 4, 5, 6}}});

      data         = {10, 20, 40, 60, 80};
      expectedData = {10 + 2, 4 * 20 + 3, 40 + 2, 4 * 60 + 3, 80 + 2};
    } else { // A Secondary rank
      data         = {20, 30, 50, 60, 70};
      expectedData = {4 * 20 + 3, 30 + 1, 50 + 2, 4 * 60 + 3, 70 + 1};
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {
      mesh->setGlobalNumberOfVertices(10);

      mesh->setVertexDistribution({{0, {1, 2, 5, 6}}, {1, {0, 1, 3, 4, 5, 7}}});

      data.assign(4, -1);
      expectedData = {2 * 20, 30, 2 * 60, 70};

    } else {
      data.assign(6, -1);
      expectedData = {10, 2 * 20, 40, 50, 2 * 60, 80};
    }
  }

  if (context.isNamed("A")) {
    c.requestConnection("B", "A");

    c.send(data);
    c.receive(data);

    BOOST_TEST(testing::equals(data, expectedData));
  } else {
    c.acceptConnection("B", "A");

    c.receive(data);
    BOOST_TEST(testing::equals(data, expectedData));
    process(data);
    c.send(data);
  }
}

/// a very similar test, but with a vertex that has been completely filtered out
void runP2PComTest2(const TestContext &context, com::PtrCommunicationFactory cf)
{
  BOOST_TEST(context.hasSize(2));
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, testing::nextMeshID()));

  m2n::PointToPointCommunication c(cf, mesh);

  vector<double> data;
  vector<double> expectedData;

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      mesh->setGlobalNumberOfVertices(10);

      mesh->setVertexDistribution({{0, {0, 1, 3, 5, 7}}, {1, {1, 2, 4, 5, 6}}});

      data         = {10, 20, 40, 60, 80};
      expectedData = {10 + 2, 4 * 20 + 3, 2 * 40 + 3, 4 * 60 + 3, 80 + 2};
    } else {
      data         = {20, 30, 50, 60, 70};
      expectedData = {4 * 20 + 3, 0, 50 + 2, 4 * 60 + 3, 70 + 1};
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {
      mesh->setGlobalNumberOfVertices(10);

      mesh->setVertexDistribution({{0, {1, 3, 5, 6}}, {1, {0, 1, 3, 4, 5, 7}}});

      data.assign(4, -1);
      expectedData = {2 * 20, 40, 2 * 60, 70};

    } else {
      data.assign(6, -1);
      expectedData = {10, 2 * 20, 40, 50, 2 * 60, 80};
    }
  }

  if (context.isNamed("A")) {
    c.requestConnection("B", "A");

    c.send(data);
    c.receive(data);
    BOOST_TEST(testing::equals(data, expectedData));
  } else {
    c.acceptConnection("B", "A");

    c.receive(data);
    BOOST_TEST(testing::equals(data, expectedData));
    process(data);
    c.send(data);
  }
}

void runSameConnectionTest(const TestContext &context, com::PtrCommunicationFactory cf)
{

  BOOST_TEST(context.hasSize(2));

  int           dimensions = 2;
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", dimensions, testing::nextMeshID()));

  if (context.isNamed("A")) {
    if (context.isPrimary()) {

      mesh->setConnectedRanks({0});
    } else {

      mesh->setConnectedRanks({1});
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {

      mesh->setConnectedRanks({0});
    } else {

      mesh->setConnectedRanks({1});
    }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  std::vector<int> receiveData;

  if (context.isNamed("A")) {
    if (context.isPrimary()) {

      c.requestPreConnection("Solid", "Fluid");
      int sendData = 5;
      c.broadcastSend(sendData);

    } else {

      c.requestPreConnection("Solid", "Fluid");
      int sendData = 10;
      c.broadcastSend(sendData);
    }
  } else {
    c.acceptPreConnection("Solid", "Fluid");
    c.broadcastReceiveAll(receiveData);

    if (context.isPrimary()) {
      BOOST_TEST(receiveData.at(0) == 5);
    } else {
      BOOST_TEST(receiveData.at(0) == 10);
    }
  }
}

void runCrossConnectionTest(const TestContext &context, com::PtrCommunicationFactory cf)
{

  BOOST_TEST(context.hasSize(2));

  int           dimensions = 2;
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", dimensions, testing::nextMeshID()));

  if (context.isNamed("A")) {
    if (context.isPrimary()) {

      mesh->setConnectedRanks({1});
    } else {

      mesh->setConnectedRanks({0});
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {

      mesh->setConnectedRanks({1});
    } else {

      mesh->setConnectedRanks({0});
    }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  std::vector<int> receiveData;

  if (context.isNamed("A")) {
    if (context.isPrimary()) {

      c.requestPreConnection("Solid", "Fluid");
      int sendData = 5;
      c.broadcastSend(sendData);

    } else {

      c.requestPreConnection("Solid", "Fluid");
      int sendData = 10;
      c.broadcastSend(sendData);
    }
  } else {
    c.acceptPreConnection("Solid", "Fluid");
    c.broadcastReceiveAll(receiveData);

    if (context.isPrimary()) {
      BOOST_TEST(receiveData.at(0) == 10);
    } else {
      BOOST_TEST(receiveData.at(0) == 5);
    }
  }
}

void runEmptyConnectionTest(const TestContext &context, com::PtrCommunicationFactory cf)
{
  BOOST_TEST(context.hasSize(2));

  int           dimensions = 2;
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", dimensions, testing::nextMeshID()));

  if (context.isNamed("A")) {
    if (context.isPrimary()) {

      mesh->setConnectedRanks({0});

    } else {
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {
      mesh->setConnectedRanks({0});

    } else {
    }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  std::vector<int> receiveData;

  if (context.isNamed("A")) {
    c.requestPreConnection("Solid", "Fluid");
    int sendData = 5;
    c.broadcastSend(sendData);
  } else {
    c.acceptPreConnection("Solid", "Fluid");
    c.broadcastReceiveAll(receiveData);

    if (context.isPrimary()) {
      BOOST_TEST(receiveData.at(0) == 5);

    } else {
      BOOST_TEST(receiveData.size() == 0);
    }
  }
}

void runP2PMeshBroadcastTest(const TestContext &context, com::PtrCommunicationFactory cf)
{
  BOOST_TEST(context.hasSize(2));

  int           dimensions = 2;
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", dimensions, testing::nextMeshID()));

  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      Eigen::VectorXd position(dimensions);
      position << 5.5, 0.0;
      mesh::Vertex &v1 = mesh->createVertex(position);
      position << 1.0, 2.0;
      mesh::Vertex &v2 = mesh->createVertex(position);
      mesh->createEdge(v1, v2);

      mesh->setConnectedRanks({0});

    } else {
      Eigen::VectorXd position(dimensions);
      position << 1.5, 0.0;
      mesh::Vertex &v1 = mesh->createVertex(position);
      position << 1.5, 2.0;
      mesh::Vertex &v2 = mesh->createVertex(position);
      mesh->createEdge(v1, v2);

      mesh->setConnectedRanks({1});
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {
      mesh->setConnectedRanks({0});

    } else {
      mesh->setConnectedRanks({1});
    }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  if (context.isNamed("A")) {

    c.requestPreConnection("Solid", "Fluid");
    c.broadcastSendMesh();
  } else {

    c.acceptPreConnection("Solid", "Fluid");
    c.broadcastReceiveAllMesh();

    if (context.isPrimary()) {
      // This rank should receive the mesh from rank 0 (fluid primary)
      BOOST_TEST(mesh->nVertices() == 2);
      BOOST_TEST(mesh->vertex(0).coord(0) == 5.50);
      BOOST_TEST(mesh->vertex(0).coord(1) == 0.0);
      BOOST_TEST(mesh->vertex(1).coord(0) == 1.0);
      BOOST_TEST(mesh->vertex(1).coord(1) == 2.0);
    } else {
      // This rank should receive the mesh from rank 1 (fluid secondary)
      BOOST_TEST(mesh->nVertices() == 2);
      BOOST_TEST(mesh->vertex(0).coord(0) == 1.50);
      BOOST_TEST(mesh->vertex(0).coord(1) == 0.0);
      BOOST_TEST(mesh->vertex(1).coord(0) == 1.50);
      BOOST_TEST(mesh->vertex(1).coord(1) == 2.0);
    }
  }
}

void runP2PComLocalCommunicationMapTest(const TestContext &context, com::PtrCommunicationFactory cf)
{
  BOOST_TEST(context.hasSize(2));

  int                             dimensions = 2;
  mesh::PtrMesh                   mesh(new mesh::Mesh("Mesh", dimensions, testing::nextMeshID()));
  const auto                      expectedId = mesh->getID();
  std::map<int, std::vector<int>> localCommunicationMap;

  if (context.isNamed("A")) {
    if (context.isPrimary()) {

      // The numbers are chosen in this way to make it easy to test weather
      // correct values are communicated or not!
      mesh->setConnectedRanks({0});
      localCommunicationMap[0].push_back(102);
      localCommunicationMap[0].push_back(1022);
      localCommunicationMap[0].push_back(10222);
      localCommunicationMap[1].push_back(103);
      localCommunicationMap[1].push_back(1033);
      localCommunicationMap[1].push_back(10333);

    } else {

      // The numbers are chosen in this way to make it easy to test weather
      // correct values are communicated or not!
      mesh->setConnectedRanks({1});
      localCommunicationMap[0].push_back(112);
      localCommunicationMap[0].push_back(1122);
      localCommunicationMap[0].push_back(11222);
      localCommunicationMap[1].push_back(113);
      localCommunicationMap[1].push_back(1133);
      localCommunicationMap[1].push_back(11333);
    }
  } else {
    BOOST_TEST(context.isNamed("B"));
    if (context.isPrimary()) {

      mesh->setConnectedRanks({0});

    } else {

      mesh->setConnectedRanks({1});
    }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  if (context.isNamed("A")) {

    c.requestPreConnection("Solid", "Fluid");
    c.scatterAllCommunicationMap(localCommunicationMap);
    BOOST_TEST(mesh->getID() == expectedId);

  } else {
    c.acceptPreConnection("Solid", "Fluid");
    c.gatherAllCommunicationMap(localCommunicationMap);
    BOOST_TEST(mesh->getID() == expectedId);
  }

  if (context.isNamed("B")) {
    if (context.isPrimary()) {
      // The numbers are chosen in this way to make it easy to test weather
      // correct values are communicated or not!
      BOOST_TEST(localCommunicationMap.size() == 1);
      BOOST_TEST(localCommunicationMap.at(0).size() == 3);
      BOOST_TEST(localCommunicationMap.at(0).at(0) == 102);
      BOOST_TEST(localCommunicationMap.at(0).at(1) == 1022);
      BOOST_TEST(localCommunicationMap.at(0).at(2) == 10222);

    } else {
      // The numbers are chosen in this way to make it easy to test weather
      // correct values are communicated or not!
      BOOST_TEST(localCommunicationMap.size() == 1);
      BOOST_TEST(localCommunicationMap.at(1).size() == 3);
      BOOST_TEST(localCommunicationMap.at(1).at(0) == 113);
      BOOST_TEST(localCommunicationMap.at(1).at(1) == 1133);
      BOOST_TEST(localCommunicationMap.at(1).at(2) == 11333);
    }
  }
}

BOOST_AUTO_TEST_SUITE(Sockets)

BOOST_AUTO_TEST_CASE(P2PComTest1)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  runP2PComTest1(context, cf);
}

BOOST_AUTO_TEST_CASE(P2PComTest2)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  runP2PComTest2(context, cf);
}

BOOST_AUTO_TEST_CASE(TestSameConnection)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  runSameConnectionTest(context, cf);
}

BOOST_AUTO_TEST_CASE(TestCrossConnection)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  runCrossConnectionTest(context, cf);
}

BOOST_AUTO_TEST_CASE(EmptyConnectionTest)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  runEmptyConnectionTest(context, cf);
}

BOOST_AUTO_TEST_CASE(P2PMeshBroadcastTest)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  runP2PMeshBroadcastTest(context, cf);
}

BOOST_AUTO_TEST_CASE(P2PComLocalCommunicationMapTest)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  runP2PComLocalCommunicationMapTest(context, cf);
}

BOOST_AUTO_TEST_SUITE_END() // Sockets

BOOST_AUTO_TEST_SUITE(MPIPorts, *boost::unit_test::label("MPI_Ports"))

BOOST_AUTO_TEST_CASE(P2PComTest1)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  runP2PComTest1(context, cf);
}

BOOST_AUTO_TEST_CASE(P2PComTest2)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  runP2PComTest2(context, cf);
}

BOOST_AUTO_TEST_CASE(TestSameConnection)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  runSameConnectionTest(context, cf);
}

BOOST_AUTO_TEST_CASE(TestCrossConnection)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  runCrossConnectionTest(context, cf);
}

BOOST_AUTO_TEST_CASE(EmptyConnectionTest)
{
  PRECICE_TEST("A"_on(2_ranks).setupIntraComm(), "B"_on(2_ranks).setupIntraComm(), Require::Events);
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  runEmptyConnectionTest(context, cf);
}

BOOST_AUTO_TEST_SUITE_END() // MPIPorts

BOOST_AUTO_TEST_SUITE_END()

#endif // not PRECICE_NO_MPI
