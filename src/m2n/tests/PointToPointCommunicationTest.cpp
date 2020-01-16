#ifndef PRECICE_NO_MPI

#include <vector>
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "m2n/PointToPointCommunication.hpp"
#include "mesh/Mesh.hpp"
#include "testing/Testing.hpp"
#include "utils/MasterSlave.hpp"

using precice::utils::MasterSlave;
using precice::utils::Parallel;

using std::vector;

using namespace precice;
using namespace m2n;

BOOST_AUTO_TEST_SUITE(M2NTests)

void process(vector<double> &data)
{
  for (auto &elem : data) {
    elem += MasterSlave::getRank() + 1;
  }
}

void P2PComTest1(com::PtrCommunicationFactory cf)
{
  BOOST_TEST(Parallel::getCommunicatorSize() == 4);

  MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true, testing::nextMeshID()));

  m2n::PointToPointCommunication c(cf, mesh);

  vector<double> data;
  vector<double> expectedData;

  switch (Parallel::getProcessRank()) {
  case 0: {
    Parallel::splitCommunicator("A.Master");

    utils::MasterSlave::configure(0, 2);

    MasterSlave::_communication->acceptConnection("A.Master", "A.Slave", "Test", 0);
    MasterSlave::_communication->setRankOffset(1);

    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(0);
    mesh->getVertexDistribution()[0].push_back(1);
    mesh->getVertexDistribution()[0].push_back(3);
    mesh->getVertexDistribution()[0].push_back(5);
    mesh->getVertexDistribution()[0].push_back(7);

    mesh->getVertexDistribution()[1].push_back(1);
    mesh->getVertexDistribution()[1].push_back(2);
    mesh->getVertexDistribution()[1].push_back(4);
    mesh->getVertexDistribution()[1].push_back(5);
    mesh->getVertexDistribution()[1].push_back(6);

    data         = {10, 20, 40, 60, 80};
    expectedData = {10 + 2, 4 * 20 + 3, 40 + 2, 4 * 60 + 3, 80 + 2};

    break;
  }
  case 1: {
    Parallel::splitCommunicator("A.Slave");
    MasterSlave::configure(1, 2);

    MasterSlave::_communication->requestConnection("A.Master", "A.Slave", "Test", 1, 1);

    data         = {20, 30, 50, 60, 70};
    expectedData = {4 * 20 + 3, 30 + 1, 50 + 2, 4 * 60 + 3, 70 + 1};

    break;
  }
  case 2: {
    Parallel::splitCommunicator("B.Master");
    MasterSlave::configure(0, 2);

    MasterSlave::_communication->acceptConnection("B.Master", "B.Slave", "Test", 0);
    MasterSlave::_communication->setRankOffset(1);

    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(1);
    mesh->getVertexDistribution()[0].push_back(2);
    mesh->getVertexDistribution()[0].push_back(5);
    mesh->getVertexDistribution()[0].push_back(6);

    mesh->getVertexDistribution()[1].push_back(0);
    mesh->getVertexDistribution()[1].push_back(1);
    mesh->getVertexDistribution()[1].push_back(3);
    mesh->getVertexDistribution()[1].push_back(4);
    mesh->getVertexDistribution()[1].push_back(5);
    mesh->getVertexDistribution()[1].push_back(7);

    data.assign(4, -1);
    expectedData = {2 * 20, 30, 2 * 60, 70};

    break;
  }
  case 3: {
    Parallel::splitCommunicator("B.Slave");
    MasterSlave::configure(1, 2);

    MasterSlave::_communication->requestConnection("B.Master", "B.Slave", "Test", 1, 1);

    data.assign(6, -1);
    expectedData = {10, 2 * 20, 40, 50, 2 * 60, 80};

    break;
  }
  }

  if (Parallel::getProcessRank() < 2) {
    c.requestConnection("B", "A");

    c.send(data.data(), data.size());
    c.receive(data.data(), data.size());

    BOOST_TEST(data == expectedData);
  } else {
    c.acceptConnection("B", "A");

    c.receive(data.data(), data.size());
    BOOST_TEST(data == expectedData);
    process(data);
    c.send(data.data(), data.size());
  }

  MasterSlave::_communication.reset();
  MasterSlave::reset();

  Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
}

/// a very similar test, but with a vertex that has been completely filtered out
void P2PComTest2(com::PtrCommunicationFactory cf)
{
  BOOST_TEST(Parallel::getCommunicatorSize() == 4);

  MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true, testing::nextMeshID()));

  m2n::PointToPointCommunication c(cf, mesh);

  vector<double> data;
  vector<double> expectedData;

  switch (Parallel::getProcessRank()) {
  case 0: {
    Parallel::splitCommunicator("A.Master");
    MasterSlave::configure(0, 2);

    MasterSlave::_communication->acceptConnection("A.Master", "A.Slave", "Test", utils::Parallel::getProcessRank());
    MasterSlave::_communication->setRankOffset(1);

    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(0);
    mesh->getVertexDistribution()[0].push_back(1);
    mesh->getVertexDistribution()[0].push_back(3);
    mesh->getVertexDistribution()[0].push_back(5);
    mesh->getVertexDistribution()[0].push_back(7);

    mesh->getVertexDistribution()[1].push_back(1);
    mesh->getVertexDistribution()[1].push_back(2);
    mesh->getVertexDistribution()[1].push_back(4);
    mesh->getVertexDistribution()[1].push_back(5);
    mesh->getVertexDistribution()[1].push_back(6);

    data         = {10, 20, 40, 60, 80};
    expectedData = {10 + 2, 4 * 20 + 3, 2 * 40 + 3, 4 * 60 + 3, 80 + 2};

    break;
  }
  case 1: {
    Parallel::splitCommunicator("A.Slave");
    utils::MasterSlave::configure(1, 2);

    MasterSlave::_communication->requestConnection("A.Master", "A.Slave", "Test", 0, 1);

    data         = {20, 30, 50, 60, 70};
    expectedData = {4 * 20 + 3, 0, 50 + 2, 4 * 60 + 3, 70 + 1};

    break;
  }
  case 2: {
    Parallel::splitCommunicator("B.Master");
    utils::MasterSlave::configure(0, 2);

    MasterSlave::_communication->acceptConnection("B.Master", "B.Slave", "Test", utils::Parallel::getProcessRank());
    MasterSlave::_communication->setRankOffset(1);

    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(1);
    mesh->getVertexDistribution()[0].push_back(3);
    mesh->getVertexDistribution()[0].push_back(5);
    mesh->getVertexDistribution()[0].push_back(6);

    mesh->getVertexDistribution()[1].push_back(0);
    mesh->getVertexDistribution()[1].push_back(1);
    mesh->getVertexDistribution()[1].push_back(3);
    mesh->getVertexDistribution()[1].push_back(4);
    mesh->getVertexDistribution()[1].push_back(5);
    mesh->getVertexDistribution()[1].push_back(7);

    data.assign(4, -1);
    expectedData = {2 * 20, 40, 2 * 60, 70};

    break;
  }
  case 3: {
    Parallel::splitCommunicator("B.Slave");
    MasterSlave::configure(1, 2);
    MasterSlave::_communication->requestConnection("B.Master", "B.Slave", "Test", 0, 1);

    data.assign(6, -1);
    expectedData = {10, 2 * 20, 40, 50, 2 * 60, 80};

    break;
  }
  }

  if (Parallel::getProcessRank() < 2) {
    c.requestConnection("B", "A");

    c.send(data.data(), data.size());
    c.receive(data.data(), data.size());
    BOOST_TEST(data == expectedData);
  } else {
    c.acceptConnection("B", "A");

    c.receive(data.data(), data.size());
    BOOST_TEST(data == expectedData);
    process(data);
    c.send(data.data(), data.size());
  }

  MasterSlave::_communication.reset();
  MasterSlave::reset();

  Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
}

void connectionTest(com::PtrCommunicationFactory cf)
{

  PRECICE_ASSERT(utils::Parallel::getCommunicatorSize() == 4);

  int           dimensions  = 2;
  bool          flipNormals = false;
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", dimensions, flipNormals, testing::nextMeshID()));

  std::vector<std::string> conections = {"same", "cross"};

  for (auto &connectionType : conections) {
    utils::MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

    switch (utils::Parallel::getProcessRank()) {
    case 0: {
      utils::Parallel::splitCommunicator("Fluid.Master");
      utils::MasterSlave::configure(0, 2);
      utils::MasterSlave::_communication->acceptConnection("Fluid.Master", "Fluid.Slave", "Test", utils::Parallel::getProcessRank());
      utils::MasterSlave::_communication->setRankOffset(1);

      if (connectionType == "same") {
        mesh->getConnectedRanks().push_back(0);
      } else {
        mesh->getConnectedRanks().push_back(1);
      }
      break;
    }
    case 1: {
      utils::Parallel::splitCommunicator("Fluid.Slave");
      utils::MasterSlave::configure(1, 2);
      utils::MasterSlave::_communication->requestConnection("Fluid.Master", "Fluid.Slave", "Test", 0, 1);

      if (connectionType == "same") {
        mesh->getConnectedRanks().push_back(1);
      } else {
        mesh->getConnectedRanks().push_back(0);
      }
      break;
    }
    case 2: {
      utils::Parallel::splitCommunicator("Solid.Master");
      utils::MasterSlave::configure(0, 2);
      utils::MasterSlave::_communication->acceptConnection("Solid.Master", "Solid.Slave", "Test", utils::Parallel::getProcessRank());
      utils::MasterSlave::_communication->setRankOffset(1);

      if (connectionType == "same") {
        mesh->getConnectedRanks().push_back(0);
      } else {
        mesh->getConnectedRanks().push_back(1);
      }
      break;
    }
    case 3: {
      utils::Parallel::splitCommunicator("Solid.Slave");
      utils::MasterSlave::configure(1, 2);
      utils::MasterSlave::_communication->requestConnection("Solid.Master", "Solid.Slave", "Test", 0, 1);

      if (connectionType == "same") {
        mesh->getConnectedRanks().push_back(1);
      } else {
        mesh->getConnectedRanks().push_back(0);
      }
      break;
    }
    }

    m2n::PointToPointCommunication c(cf, mesh);

    std::vector<int> receiveData;

    if (utils::Parallel::getProcessRank() == 0) {

      c.requestPreConnection("Solid", "Fluid");
      int sendData = 5;
      c.broadcastSend(sendData);

    } else if (utils::Parallel::getProcessRank() == 1) {

      c.requestPreConnection("Solid", "Fluid");
      int sendData = 10;
      c.broadcastSend(sendData);

    } else {
      c.acceptPreConnection("Solid", "Fluid");
      c.broadcastReceiveAll(receiveData);
    }

    if (utils::Parallel::getProcessRank() == 2) {
      if (connectionType == "same") {
        BOOST_TEST(receiveData[0] == 5);
      } else {
        BOOST_TEST(receiveData[1] == 10);
      }

    } else if (utils::Parallel::getProcessRank() == 3) {
      if (connectionType == "same") {
        BOOST_TEST(receiveData[0] == 10);
      } else {
        BOOST_TEST(receiveData[1] == 5);
      }
    }

    utils::MasterSlave::_communication = nullptr;
    utils::MasterSlave::reset();
    utils::Parallel::synchronizeProcesses();
    utils::Parallel::clearGroups();
    mesh::Data::resetDataCount();
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
  }
}

void emptyConnectionTest(com::PtrCommunicationFactory cf)
{

  PRECICE_ASSERT(utils::Parallel::getCommunicatorSize() == 4);

  int           dimensions  = 2;
  bool          flipNormals = false;
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", dimensions, flipNormals, testing::nextMeshID()));

  utils::MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    utils::Parallel::splitCommunicator("Fluid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Fluid.Master", "Fluid.Slave", "Test", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);

    mesh->getConnectedRanks().push_back(0);

    break;
  }
  case 1: {
    utils::Parallel::splitCommunicator("Fluid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Fluid.Master", "Fluid.Slave", "Test", 0, 1);

    break;
  }
  case 2: {
    utils::Parallel::splitCommunicator("Solid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Solid.Master", "Solid.Slave", "Test", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);

    mesh->getConnectedRanks().push_back(0);

    break;
  }
  case 3: {
    utils::Parallel::splitCommunicator("Solid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Solid.Master", "Solid.Slave", "Test", 0, 1);

    break;
  }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  std::vector<int> receiveData;

  if (utils::Parallel::getProcessRank() < 2) {

    c.requestPreConnection("Solid", "Fluid");
    int sendData = 5;
    c.broadcastSend(sendData);

  } else if (utils::Parallel::getProcessRank() > 1) {
    c.acceptPreConnection("Solid", "Fluid");
    c.broadcastReceiveAll(receiveData);
  }

  if (utils::Parallel::getProcessRank() == 2) {
    BOOST_TEST(receiveData[0] == 5);

  } else if (utils::Parallel::getProcessRank() == 3) {
    BOOST_TEST(receiveData.size() == 0);
  }

  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
}

void P2PMeshBroadcastTest(com::PtrCommunicationFactory cf)
{
  PRECICE_ASSERT(utils::Parallel::getCommunicatorSize() == 4);
  utils::MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  int           dimensions  = 2;
  bool          flipNormals = false;
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", dimensions, flipNormals, testing::nextMeshID()));

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    utils::Parallel::splitCommunicator("Fluid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Fluid.Master", "Fluid.Slave", "Test", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);

    Eigen::VectorXd position(dimensions);
    position << 5.5, 0.0;
    mesh::Vertex &v1 = mesh->createVertex(position);
    position << 1.0, 2.0;
    mesh::Vertex &v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);

    mesh->getConnectedRanks().push_back(0);

    break;
  }
  case 1: {
    utils::Parallel::splitCommunicator("Fluid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Fluid.Master", "Fluid.Slave", "Test", 0, 1);

    Eigen::VectorXd position(dimensions);
    position << 1.5, 0.0;
    mesh::Vertex &v1 = mesh->createVertex(position);
    position << 1.5, 2.0;
    mesh::Vertex &v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);

    mesh->getConnectedRanks().push_back(1);

    break;
  }
  case 2: {
    utils::Parallel::splitCommunicator("Solid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Solid.Master", "Solid.Slave", "Test", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);

    mesh->getConnectedRanks().push_back(0);

    break;
  }
  case 3: {
    utils::Parallel::splitCommunicator("Solid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Solid.Master", "Solid.Slave", "Test", 0, 1);

    mesh->getConnectedRanks().push_back(1);

    break;
  }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  if (utils::Parallel::getProcessRank() < 2) {

    c.requestPreConnection("Solid", "Fluid");
    c.broadcastSendMesh();
  } else {

    c.acceptPreConnection("Solid", "Fluid");
    c.broadcastReceiveMesh();

    if (utils::Parallel::getProcessRank() == 2) {
      // This rank should receive the mesh from rank 0 (fluid master)
      BOOST_TEST(mesh->vertices().size() == 2);
      BOOST_TEST(mesh->vertices()[0].getCoords()[0] == 5.50);
      BOOST_TEST(mesh->vertices()[0].getCoords()[1] == 0.0);
      BOOST_TEST(mesh->vertices()[1].getCoords()[0] == 1.0);
      BOOST_TEST(mesh->vertices()[1].getCoords()[1] == 2.0);
    }

    if (utils::Parallel::getProcessRank() == 3) {
      // This rank should receive the mesh from rank 1 (fluid slave)
      BOOST_TEST(mesh->vertices().size() == 2);
      BOOST_TEST(mesh->vertices()[0].getCoords()[0] == 1.50);
      BOOST_TEST(mesh->vertices()[0].getCoords()[1] == 0.0);
      BOOST_TEST(mesh->vertices()[1].getCoords()[0] == 1.50);
      BOOST_TEST(mesh->vertices()[1].getCoords()[1] == 2.0);
    }
  }

  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
}

void P2PComLCMTest(com::PtrCommunicationFactory cf)
{
  PRECICE_ASSERT(utils::Parallel::getCommunicatorSize() == 4);
  utils::MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  int                             dimensions  = 2;
  bool                            flipNormals = false;
  mesh::PtrMesh                   mesh(new mesh::Mesh("Mesh", dimensions, flipNormals, testing::nextMeshID()));
  const auto                      expectedId = mesh->getID();
  std::map<int, std::vector<int>> localCommunicationMap;

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    utils::Parallel::splitCommunicator("Fluid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Fluid.Master", "Fluid.Slave", "Test", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);

    // The numbers are chosen in this way to make it easy to test weather
    // correct values are communicated or not!
    mesh->getConnectedRanks().push_back(0);
    localCommunicationMap[0].push_back(102);
    localCommunicationMap[0].push_back(1022);
    localCommunicationMap[0].push_back(10222);
    localCommunicationMap[1].push_back(103);
    localCommunicationMap[1].push_back(1033);
    localCommunicationMap[1].push_back(10333);

    break;
  }
  case 1: {
    utils::Parallel::splitCommunicator("Fluid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Fluid.Master", "Fluid.Slave", "Test", 0, 1);

    // The numbers are chosen in this way to make it easy to test weather
    // correct values are communicated or not!
    mesh->getConnectedRanks().push_back(1);
    localCommunicationMap[0].push_back(112);
    localCommunicationMap[0].push_back(1122);
    localCommunicationMap[0].push_back(11222);
    localCommunicationMap[1].push_back(113);
    localCommunicationMap[1].push_back(1133);
    localCommunicationMap[1].push_back(11333);

    break;
  }
  case 2: {
    utils::Parallel::splitCommunicator("Solid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Solid.Master", "Solid.Slave", "Test", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);

    mesh->getConnectedRanks().push_back(0);

    break;
  }
  case 3: {
    utils::Parallel::splitCommunicator("Solid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Solid.Master", "Solid.Slave", "Test", 0, 1);

    mesh->getConnectedRanks().push_back(1);

    break;
  }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  if (utils::Parallel::getProcessRank() < 2) {

    c.requestPreConnection("Solid", "Fluid");
    c.broadcastSendLCM(localCommunicationMap);
    BOOST_TEST(mesh->getID() == expectedId);

  } else {
    c.acceptPreConnection("Solid", "Fluid");
    c.broadcastReceiveLCM(localCommunicationMap);
    BOOST_TEST(mesh->getID() == expectedId);
  }

  if (utils::Parallel::getProcessRank() == 2) {
    // The numbers are chosen in this way to make it easy to test weather
    // correct values are communicated or not!
    BOOST_TEST(localCommunicationMap.size() == 1);
    BOOST_TEST(localCommunicationMap[0].size() == 3);
    BOOST_TEST(localCommunicationMap[0][0] == 102);
    BOOST_TEST(localCommunicationMap[0][1] == 1022);
    BOOST_TEST(localCommunicationMap[0][2] == 10222);

  } else if (utils::Parallel::getProcessRank() == 3) {
    // The numbers are chosen in this way to make it easy to test weather
    // correct values are communicated or not!
    BOOST_TEST(localCommunicationMap.size() == 1);
    BOOST_TEST(localCommunicationMap[1].size() == 3);
    BOOST_TEST(localCommunicationMap[1][0] == 113);
    BOOST_TEST(localCommunicationMap[1][1] == 1133);
    BOOST_TEST(localCommunicationMap[1][2] == 11333);
  }

  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
}

BOOST_AUTO_TEST_CASE(SocketCommunication,
                     *testing::OnSize(4))
{
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  if (utils::Parallel::getProcessRank() < 4) {
    P2PComTest1(cf);
    P2PComTest2(cf);
    connectionTest(cf);
    emptyConnectionTest(cf);
    P2PMeshBroadcastTest(cf);
    P2PComLCMTest(cf);
  }
}

BOOST_AUTO_TEST_CASE(MPIPortsCommunication,
                     *testing::OnSize(4) * boost::unit_test::label("MPI_Ports"))
{
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  if (utils::Parallel::getProcessRank() < 4) {
    P2PComTest1(cf);
    P2PComTest2(cf);
    connectionTest(cf);
    emptyConnectionTest(cf);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif // not PRECICE_NO_MPI
