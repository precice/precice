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

    MasterSlave::_rank       = 1;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = false;
    MasterSlave::_slaveMode  = true;

    MasterSlave::_communication->requestConnection("A.Master", "A.Slave", 1, 1);

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

    MasterSlave::_rank       = 1;
    MasterSlave::_size       = 2;
    MasterSlave::_masterMode = false;
    MasterSlave::_slaveMode  = true;

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

BOOST_AUTO_TEST_CASE(SocketCommunication,
                     * testing::OnSize(4))
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

BOOST_AUTO_TEST_SUITE_END()

#endif // not PRECICE_NO_MPI
