// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#include "PointToPointCommunicationTest.hpp"

#include "m2n/PointToPointCommunication.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "mesh/Mesh.hpp"

#include "tarch/tests/TestCaseFactory.h"

#include <vector>

using precice::utils::Parallel;
using precice::utils::MasterSlave;

using std::vector;
using std::rand;

registerTest(precice::m2n::tests::PointToPointCommunicationTest);

namespace precice {
namespace m2n {
namespace tests {

void
process(vector<double>& data) {
  for (int i = 0; i < data.size(); ++i) {
    data[i] += MasterSlave::_rank + 1;
  }
}

bool
equal(vector<double> const& data, vector<double> const& expectedData) {
  bool valid = true;

  if (data.size() != expectedData.size())
    return false;

  for (int i = 0; i < data.size(); ++i) {
    valid &= (data[i] == expectedData[i]);
  }

  return valid;
}

tarch::logging::Log PointToPointCommunicationTest::_log(
    "precice::m2n::tests::PointToPointCommunicationTest");

PointToPointCommunicationTest::PointToPointCommunicationTest()
    : TestCase("m2n::tests::PointToPointCommunicationTest") {
}

void
PointToPointCommunicationTest::run() {
  preciceTrace("run");

  Parallel::synchronizeProcesses();

  if (Parallel::getCommunicatorSize() > 3) {
    auto communicator = Parallel::getRestrictedCommunicator({0, 1, 2, 3});

    if (Parallel::getProcessRank() < 4) {
      Parallel::setGlobalCommunicator(communicator);
      testMethod(testSocketCommunication);
      testMethod(testMPIPortsCommunication);
      Parallel::setGlobalCommunicator(Parallel::getCommunicatorWorld());
    }
  }
}

void
PointToPointCommunicationTest::testSocketCommunication() {
  preciceTrace("testSocketCommunication");

  com::CommunicationFactory::SharedPointer cf(
      new com::SocketCommunicationFactory);

  test(cf);
}

void
PointToPointCommunicationTest::testMPIPortsCommunication() {
  preciceTrace("testMPIDirectCommunication");

  com::CommunicationFactory::SharedPointer cf(
      new com::MPIPortsCommunicationFactory);

  test(cf);
}

void
PointToPointCommunicationTest::test(
    com::CommunicationFactory::SharedPointer cf) {
  assertion(Parallel::getCommunicatorSize() == 4);

  validateEquals(Parallel::getCommunicatorSize(), 4);

  MasterSlave::_communication =
      com::Communication::SharedPointer(new com::MPIDirectCommunication);

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));

  m2n::PointToPointCommunication c(cf, mesh);

  vector<double> data;
  vector<double> expectedData;

  switch (Parallel::getProcessRank()) {
  case 0: {
    Parallel::initialize(NULL, NULL, "A.Master");

    MasterSlave::_rank = 0;
    MasterSlave::_size = 2;
    MasterSlave::_masterMode = true;
    MasterSlave::_slaveMode = false;

    MasterSlave::_communication->acceptConnection("A.Master", "A.Slave", 0, 1);
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

    data = {10, 20, 40, 60, 80};
    expectedData = {10 + 2, 4 * 20 + 3, 40 + 2, 4 * 60 + 3, 80 + 2};

    break;
  }
  case 1: {
    Parallel::initialize(NULL, NULL, "A.Slave");

    MasterSlave::_rank = 1;
    MasterSlave::_size = 2;
    MasterSlave::_masterMode = false;
    MasterSlave::_slaveMode = true;

    MasterSlave::_communication->requestConnection("A.Master", "A.Slave", 0, 1);

    data = {20, 30, 50, 60, 70};
    expectedData = {4 * 20 + 3, 30 + 1, 50 + 2, 4 * 60 + 3, 70 + 1};

    break;
  }
  case 2: {
    Parallel::initialize(NULL, NULL, "B.Master");

    MasterSlave::_rank = 0;
    MasterSlave::_size = 2;
    MasterSlave::_masterMode = true;
    MasterSlave::_slaveMode = false;

    MasterSlave::_communication->acceptConnection("B.Master", "B.Slave", 0, 1);
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

    data = {rand(), rand(), rand(), rand()};
    expectedData = {2 * 20, 30, 2 * 60, 70};

    break;
  }
  case 3: {
    Parallel::initialize(NULL, NULL, "B.Slave");

    MasterSlave::_rank = 1;
    MasterSlave::_size = 2;
    MasterSlave::_masterMode = false;
    MasterSlave::_slaveMode = true;

    MasterSlave::_communication->requestConnection("B.Master", "B.Slave", 0, 1);

    data = {rand(), rand(), rand(), rand(), rand(), rand()};
    expectedData = {10, 2 * 20, 40, 50, 2 * 60, 80};

    break;
  }
  }

  if (Parallel::getProcessRank() < 2) {
    c.requestConnection("B", "A");

    c.send(data.data(), data.size());

    c.receive(data.data(), data.size());

    validate(equal(data, expectedData));
  } else {
    c.acceptConnection("B", "A");

    c.receive(data.data(), data.size());

    validate(equal(data, expectedData));

    process(data);

    c.send(data.data(), data.size());
  }

  MasterSlave::_communication.reset();
  MasterSlave::_rank = Parallel::getProcessRank();
  MasterSlave::_size = Parallel::getCommunicatorSize();
  MasterSlave::_masterMode = false;
  MasterSlave::_slaveMode = false;

  Parallel::synchronizeProcesses();
}
}
}
} // namespace precice, m2n, tests

#endif // not PRECICE_NO_MPI
