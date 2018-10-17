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
    elem += MasterSlave::_rank + 1;
  }
}

void P2PComTest1(com::PtrCommunicationFactory cf)
{
  assertion(Parallel::getCommunicatorSize() == 4);

  MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));

  m2n::PointToPointCommunication c(cf, mesh);

  vector<double> data;
  vector<double> expectedData;

  switch (Parallel::getProcessRank()) {
  case 0: {
    Parallel::splitCommunicator("A.Master");

    MasterSlave::_rank       = 0;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = true;
    MasterSlave::_slaveMode  = false;

    MasterSlave::_communication->acceptConnection("A.Master", "A.Slave");
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

    MasterSlave::_rank       = 1;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = false;
    MasterSlave::_slaveMode  = true;

    MasterSlave::_communication->requestConnection("A.Master", "A.Slave", 0, 1);

    data         = {20, 30, 50, 60, 70};
    expectedData = {4 * 20 + 3, 30 + 1, 50 + 2, 4 * 60 + 3, 70 + 1};

    break;
  }
  case 2: {
    Parallel::splitCommunicator("B.Master");

    MasterSlave::_rank       = 0;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = true;
    MasterSlave::_slaveMode  = false;

    MasterSlave::_communication->acceptConnection("B.Master", "B.Slave");
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

    MasterSlave::_rank       = 1;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = false;
    MasterSlave::_slaveMode  = true;

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

/// a very similar test, but with a vertex that has been completely filtered out
void P2PComTest2(com::PtrCommunicationFactory cf)
{
  assertion(Parallel::getCommunicatorSize() == 4);

  MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));

  m2n::PointToPointCommunication c(cf, mesh);

  vector<double> data;
  vector<double> expectedData;

  switch (Parallel::getProcessRank()) {
  case 0: {
    Parallel::splitCommunicator("A.Master");

    MasterSlave::_rank       = 0;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = true;
    MasterSlave::_slaveMode  = false;

    MasterSlave::_communication->acceptConnection("A.Master", "A.Slave");
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

    MasterSlave::_rank       = 1;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = false;
    MasterSlave::_slaveMode  = true;

    MasterSlave::_communication->requestConnection("A.Master", "A.Slave", 0, 1);

    data         = {20, 30, 50, 60, 70};
    expectedData = {4 * 20 + 3, 0, 50 + 2, 4 * 60 + 3, 70 + 1};

    break;
  }
  case 2: {
    Parallel::splitCommunicator("B.Master");

    MasterSlave::_rank       = 0;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = true;
    MasterSlave::_slaveMode  = false;

    MasterSlave::_communication->acceptConnection("B.Master", "B.Slave");
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

    MasterSlave::_rank       = 1;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = false;
    MasterSlave::_slaveMode  = true;

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

BOOST_AUTO_TEST_CASE(SocketCommunication, * testing::OnSize(4))
{
  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  if (utils::Parallel::getProcessRank() < 4) {
    P2PComTest1(cf);
    P2PComTest2(cf);
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
  }
}

BOOST_AUTO_TEST_CASE(P2PComMeshTest, * testing::OnSize(4))
{
  
  assertion(utils::Parallel::getCommunicatorSize() == 4);

//  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
//    com::PtrCommunicationFactory cf =  com::PtrCommunicationFactory(new com::SocketCommunicationFactory);
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  
  utils::MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));  

  int dimensions = 2;
  bool flipNormals = true;


  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    utils::Parallel::splitCommunicator("Fluid.Master");

    utils::MasterSlave::_rank       = 0;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = true;
    utils::MasterSlave::_slaveMode  = false;

    utils::MasterSlave::_communication->acceptConnection("Fluid.Master", "Fluid.Slave");
    utils::MasterSlave::_communication->setRankOffset(1);

    Eigen::VectorXd position(dimensions);
    position <<0.5, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 1.0, 2.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);

    mesh->getCommunicationMap()[0].push_back(-1);
//    mesh->getCommunicationMap()[1].push_back(-1);
    
    
    break;
  }
  case 1: {
    utils::Parallel::splitCommunicator("Fluid.Slave");

    utils::MasterSlave::_rank       = 1;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = false;
    utils::MasterSlave::_slaveMode  = true;

    utils::MasterSlave::_communication->requestConnection("Fluid.Master", "Fluid.Slave", 0, 1);

    Eigen::VectorXd position(dimensions);
    position <<1.5, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 1.5, 2.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);

    //   mesh->getCommunicationMap()[0].push_back(-1);
    mesh->getCommunicationMap()[1].push_back(-1);
    


    break;
  }
  case 2: {
    utils::Parallel::splitCommunicator("Solid.Master");

    utils::MasterSlave::_rank       = 0;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = true;
    utils::MasterSlave::_slaveMode  = false;

    utils::MasterSlave::_communication->acceptConnection("Solid.Master", "Solid.Slave");
    utils::MasterSlave::_communication->setRankOffset(1);


    Eigen::VectorXd position(dimensions);
    position << 1.5, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 1.0, 2.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);

    mesh->getCommunicationMap()[0].push_back(-1);
    //mesh->getCommunicationMap()[1].push_back(-1);
    

    break;
  }
  case 3: {
    utils::Parallel::splitCommunicator("Solid.Slave");

    utils::MasterSlave::_rank       = 1;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = false;
    utils::MasterSlave::_slaveMode  = true;

    utils::MasterSlave::_communication->requestConnection("Solid.Master", "Solid.Slave", 0, 1);

    Eigen::VectorXd position(dimensions);
    position <<1.0, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 1.0, 2.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);
    
    //mesh->getCommunicationMap()[0].push_back(-1);
    mesh->getCommunicationMap()[1].push_back(-1);


    break;
  }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  if (utils::Parallel::getProcessRank() < 2) {
  
    c.requestPreConnection("Solid", "Fluid");
    c.sendMesh(*mesh);
  } else {

    c.acceptPreConnection("Solid", "Fluid");
    c.receiveMesh(*mesh);

      if(utils::Parallel::getProcessRank() ==2 )
      {        
        BOOST_TEST(mesh->vertices()[2].getCoords()[0]==0.50);
        BOOST_TEST(mesh->vertices().size()==4);
      }

      if(utils::Parallel::getProcessRank() ==3 )
      {       
        BOOST_TEST(mesh->vertices()[2].getCoords()[0]==1.50);
        BOOST_TEST(mesh->vertices().size()==4);    
      }
    
  }
  
  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
}


BOOST_AUTO_TEST_CASE(P2PComLCMTest, * testing::OnSize(4))
{
  
  assertion(utils::Parallel::getCommunicatorSize() == 4);

//  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  
  utils::MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));
  mesh::Mesh::FeedbackMap localCommunicationMap;

  int dimensions = 2;
  bool flipNormals = true;


  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    utils::Parallel::splitCommunicator("Fluid.Master");

    utils::MasterSlave::_rank       = 0;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = true;
    utils::MasterSlave::_slaveMode  = false;

    utils::MasterSlave::_communication->acceptConnection("Fluid.Master", "Fluid.Slave");
    utils::MasterSlave::_communication->setRankOffset(1);
    

    mesh->getCommunicationMap()[0].push_back(-1);
    mesh->getCommunicationMap()[1].push_back(-1);
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

    utils::MasterSlave::_rank       = 1;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = false;
    utils::MasterSlave::_slaveMode  = true;

    utils::MasterSlave::_communication->requestConnection("Fluid.Master", "Fluid.Slave", 0, 1);
    

    mesh->getCommunicationMap()[0].push_back(-1);
    mesh->getCommunicationMap()[1].push_back(-1);
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

    utils::MasterSlave::_rank       = 0;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = true;
    utils::MasterSlave::_slaveMode  = false;

    utils::MasterSlave::_communication->acceptConnection("Solid.Master", "Solid.Slave");
    utils::MasterSlave::_communication->setRankOffset(1);
    

    mesh->getCommunicationMap()[0].push_back(-1);
    mesh->getCommunicationMap()[1].push_back(-1);
    

    break;
  }
  case 3: {
    utils::Parallel::splitCommunicator("Solid.Slave");

    utils::MasterSlave::_rank       = 1;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = false;
    utils::MasterSlave::_slaveMode  = true;

    utils::MasterSlave::_communication->requestConnection("Solid.Master", "Solid.Slave", 0, 1);


    mesh->getCommunicationMap()[0].push_back(-1);
    mesh->getCommunicationMap()[1].push_back(-1);


    break;
  }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  if (utils::Parallel::getProcessRank() < 2) {
  
    c.requestPreConnection("Solid", "Fluid");
    c.sendCommunicationMap(localCommunicationMap);
    BOOST_TEST(mesh->getID()==0);
   
  } else
  {
    c.acceptPreConnection("Solid", "Fluid");
    c.receiveCommunicationMap(localCommunicationMap);
    BOOST_TEST(mesh->getID()==0);
  }

 if(utils::Parallel::getProcessRank() == 2 )
  {    
    BOOST_TEST(localCommunicationMap.size() == 2);
    BOOST_TEST(localCommunicationMap[0].size() ==3);
    BOOST_TEST(localCommunicationMap[1].size() ==3);
    BOOST_TEST(localCommunicationMap[0][0] ==102);
    BOOST_TEST(localCommunicationMap[0][1] ==1022);
    BOOST_TEST(localCommunicationMap[0][2] ==10222);
    BOOST_TEST(localCommunicationMap[1][0] ==112);
    BOOST_TEST(localCommunicationMap[1][1] ==1122);
    BOOST_TEST(localCommunicationMap[1][2] ==11222);                   
  } else if(utils::Parallel::getProcessRank() == 3 )
  {
    BOOST_TEST(localCommunicationMap.size() == 2);
    BOOST_TEST(localCommunicationMap[0].size() ==3);
    BOOST_TEST(localCommunicationMap[1].size() ==3);
    BOOST_TEST(localCommunicationMap[0][0] ==103);
    BOOST_TEST(localCommunicationMap[0][1] ==1033);
    BOOST_TEST(localCommunicationMap[0][2] ==10333);
    BOOST_TEST(localCommunicationMap[1][0] ==113);
    BOOST_TEST(localCommunicationMap[1][1] ==1133);
    BOOST_TEST(localCommunicationMap[1][2] ==11333);   
  }
  
  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
  
}

BOOST_AUTO_TEST_SUITE_END()

#endif // not PRECICE_NO_MPI
