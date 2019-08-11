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

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));

  m2n::PointToPointCommunication c(cf, mesh);

  vector<double> data;
  vector<double> expectedData;

  switch (Parallel::getProcessRank()) {
  case 0: {
    Parallel::splitCommunicator("A.Master");

    utils::MasterSlave::configure(0, 2);

    MasterSlave::_communication->acceptConnection("A.Master", "A.Slave", 0);
    MasterSlave::_communication->setRankOffset(1);

    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(0);
    mesh->getVertexDistribution()[0].push_back(1); // <-
    mesh->getVertexDistribution()[0].push_back(3);
    mesh->getVertexDistribution()[0].push_back(5); // <-
    mesh->getVertexDistribution()[0].push_back(7);

    mesh->getVertexDistribution()[1].push_back(1); // <-
    mesh->getVertexDistribution()[1].push_back(2);
    mesh->getVertexDistribution()[1].push_back(4);
    mesh->getVertexDistribution()[1].push_back(5); // <-
    mesh->getVertexDistribution()[1].push_back(6);

    data         = {10, 20, 40, 60, 80};
    expectedData = {10 + 2, 4 * 20 + 3, 40 + 2, 4 * 60 + 3, 80 + 2};

    break;
  }
  case 1: {
    Parallel::splitCommunicator("A.Slave");
    MasterSlave::configure(1, 2);

    MasterSlave::_communication->requestConnection("A.Master", "A.Slave", 1, 1);

    data         = {20, 30, 50, 60, 70};
    expectedData = {4 * 20 + 3, 30 + 1, 50 + 2, 4 * 60 + 3, 70 + 1};

    break;
  }
  case 2: {
    Parallel::splitCommunicator("B.Master");
    MasterSlave::configure(0, 2);

    MasterSlave::_communication->acceptConnection("B.Master", "B.Slave", 0);
    MasterSlave::_communication->setRankOffset(1);

    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(1); // <-
    mesh->getVertexDistribution()[0].push_back(2);
    mesh->getVertexDistribution()[0].push_back(5); // <-
    mesh->getVertexDistribution()[0].push_back(6);

    mesh->getVertexDistribution()[1].push_back(0);
    mesh->getVertexDistribution()[1].push_back(1); // <-
    mesh->getVertexDistribution()[1].push_back(3);
    mesh->getVertexDistribution()[1].push_back(4);
    mesh->getVertexDistribution()[1].push_back(5); // <-
    mesh->getVertexDistribution()[1].push_back(7);

    data.assign(4, -1);
    expectedData = {2 * 20, 30, 2 * 60, 70};

    break;
  }
  case 3: {
    Parallel::splitCommunicator("B.Slave");
    MasterSlave::configure(1, 2);

    MasterSlave::_communication->requestConnection("B.Master", "B.Slave", 1, 1);

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

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));

  m2n::PointToPointCommunication c(cf, mesh);

  vector<double> data;
  vector<double> expectedData;

  switch (Parallel::getProcessRank()) {
  case 0: {
    Parallel::splitCommunicator("A.Master");
    MasterSlave::configure(0, 2);

    MasterSlave::_communication->acceptConnection("A.Master", "A.Slave", utils::Parallel::getProcessRank());
    MasterSlave::_communication->setRankOffset(1);

    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(0);
    mesh->getVertexDistribution()[0].push_back(1); // <-
    mesh->getVertexDistribution()[0].push_back(3);
    mesh->getVertexDistribution()[0].push_back(5); // <-
    mesh->getVertexDistribution()[0].push_back(7);

    mesh->getVertexDistribution()[1].push_back(1); // <-
    mesh->getVertexDistribution()[1].push_back(2);
    mesh->getVertexDistribution()[1].push_back(4);
    mesh->getVertexDistribution()[1].push_back(5); // <-
    mesh->getVertexDistribution()[1].push_back(6);

    data         = {10, 20, 40, 60, 80};
    expectedData = {10 + 2, 4 * 20 + 3, 2 * 40 + 3, 4 * 60 + 3, 80 + 2};

    break;
  }
  case 1: {
    Parallel::splitCommunicator("A.Slave");
    utils::MasterSlave::configure(1, 2);

    MasterSlave::_communication->requestConnection("A.Master", "A.Slave", 0, 1);

    data         = {20, 30, 50, 60, 70};
    expectedData = {4 * 20 + 3, 0, 50 + 2, 4 * 60 + 3, 70 + 1};

    break;
  }
  case 2: {
    Parallel::splitCommunicator("B.Master");
    utils::MasterSlave::configure(0, 2);

    MasterSlave::_communication->acceptConnection("B.Master", "B.Slave", utils::Parallel::getProcessRank());
    MasterSlave::_communication->setRankOffset(1);

    mesh->setGlobalNumberOfVertices(10);

    mesh->getVertexDistribution()[0].push_back(1); // <-
    mesh->getVertexDistribution()[0].push_back(3);
    mesh->getVertexDistribution()[0].push_back(5); // <-
    mesh->getVertexDistribution()[0].push_back(6);

    mesh->getVertexDistribution()[1].push_back(0);
    mesh->getVertexDistribution()[1].push_back(1); // <-
    mesh->getVertexDistribution()[1].push_back(3);
    mesh->getVertexDistribution()[1].push_back(4);
    mesh->getVertexDistribution()[1].push_back(5); // <-
    mesh->getVertexDistribution()[1].push_back(7);

    data.assign(4, -1);
    expectedData = {2 * 20, 40, 2 * 60, 70};

    break;
  }
  case 3: {
    Parallel::splitCommunicator("B.Slave");
    MasterSlave::configure(1, 2);
    MasterSlave::_communication->requestConnection("B.Master", "B.Slave", 0, 1);

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

void connectionTest1(com::PtrCommunicationFactory cf)
{

  assertion(utils::Parallel::getCommunicatorSize() == 4);  
  utils::MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();
  
  int dimensions = 2;
  bool flipNormals = false;
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", dimensions, flipNormals));


  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    
    utils::Parallel::splitCommunicator("Fluid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Fluid.Master", "Fluid.Slave", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);
  
    mesh->getConnectedRanks().push_back(0);
   
    break;
  }
  case 1: {
    utils::Parallel::splitCommunicator("Fluid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Fluid.Master", "Fluid.Slave", 0, 1);
  
    mesh->getConnectedRanks().push_back(1);
    break;
  }
  case 2: {
    utils::Parallel::splitCommunicator("Solid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Solid.Master", "Solid.Slave", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);
    
    mesh->getConnectedRanks().push_back(0);
    
    break;
  }
  case 3: {
    utils::Parallel::splitCommunicator("Solid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Solid.Master", "Solid.Slave", 0, 1);
    
    mesh->getConnectedRanks().push_back(1);

    break;
  }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  double receiveData = 0;

  if (utils::Parallel::getProcessRank() == 0) {
  
    c.requestPreConnection("Solid", "Fluid");
    double sendData = 5;
    c.broadcastSend(sendData);    
   
  } else if (utils::Parallel::getProcessRank() == 1) {
  
    c.requestPreConnection("Solid", "Fluid");
    double sendData = 10;
    c.broadcastSend(sendData);    
   
  } else
  {    
    c.acceptPreConnection("Solid", "Fluid");
    c.broadcastReceive(receiveData);    
  }

  if(utils::Parallel::getProcessRank() == 2 )
  {

    BOOST_TEST(receiveData == 5);
    
  } else if(utils::Parallel::getProcessRank() == 3 )
  {
    
    BOOST_TEST(receiveData == 10);
  }
  
  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());

}

void connectionTest2(com::PtrCommunicationFactory cf)
{

  assertion(utils::Parallel::getCommunicatorSize() == 4);  
  utils::MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  int dimensions = 2;
  bool flipNormals = false;
  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", dimensions, flipNormals));

  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    
    utils::Parallel::splitCommunicator("Fluid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Fluid.Master", "Fluid.Slave", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);
  
    mesh->getConnectedRanks().push_back(1);
   
    break;
  }
  case 1: {
    utils::Parallel::splitCommunicator("Fluid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Fluid.Master", "Fluid.Slave", 0, 1);
  
    mesh->getConnectedRanks().push_back(0);
    break;
  }
  case 2: {
    utils::Parallel::splitCommunicator("Solid.Master");
    utils::MasterSlave::configure(0, 2);
    utils::MasterSlave::_communication->acceptConnection("Solid.Master", "Solid.Slave", utils::Parallel::getProcessRank());
    utils::MasterSlave::_communication->setRankOffset(1);
    
    mesh->getConnectedRanks().push_back(1);
    
    break;
  }
  case 3: {
    utils::Parallel::splitCommunicator("Solid.Slave");
    utils::MasterSlave::configure(1, 2);
    utils::MasterSlave::_communication->requestConnection("Solid.Master", "Solid.Slave", 0, 1);
    
    mesh->getConnectedRanks().push_back(0);

    break;
  }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  double receiveData = 0;

  if (utils::Parallel::getProcessRank() == 0) {
  
    c.requestPreConnection("Solid", "Fluid");
    double sendData = 5;
    c.broadcastSend(sendData);    
   
  } else if (utils::Parallel::getProcessRank() == 1) {
  
    c.requestPreConnection("Solid", "Fluid");
    double sendData = 10;
    c.broadcastSend(sendData);    
   
  } else
  {    
    c.acceptPreConnection("Solid", "Fluid");
    c.broadcastReceive(receiveData);    
  }

  if(utils::Parallel::getProcessRank() == 2 )
  {

    BOOST_TEST(receiveData == 10);
    
  } else if(utils::Parallel::getProcessRank() == 3 )
  {
    
    BOOST_TEST(receiveData == 5);
  }
  
  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());

}

BOOST_AUTO_TEST_CASE(SocketCommunication,
                     * testing::OnSize(4))
{
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  if (utils::Parallel::getProcessRank() < 4) {
    P2PComTest1(cf);
    P2PComTest2(cf);
    connectionTest1(cf);
    connectionTest2(cf);
  }
}

BOOST_AUTO_TEST_CASE(MPIPortsCommunication,
                     * testing::OnSize(4)
                     * boost::unit_test::label("MPI_Ports"))
{
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  if (utils::Parallel::getProcessRank() < 4) {
    P2PComTest1(cf);
    P2PComTest2(cf);
    connectionTest1(cf);
    connectionTest2(cf);
  }
}


BOOST_AUTO_TEST_SUITE_END()

#endif // not PRECICE_NO_MPI
