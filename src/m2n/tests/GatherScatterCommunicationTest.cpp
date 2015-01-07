// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NO_MPI

#include "GatherScatterCommunicationTest.hpp"
#include "utils/Parallel.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/GatherScatterCommunication.hpp"
#include "m2n/SharedPointer.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include "utils/MasterSlave.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::m2n::tests::GatherScatterCommunicationTest)

namespace precice {
namespace m2n {
namespace tests {

tarch::logging::Log GatherScatterCommunicationTest::
  _log ( "precice::m2n::tests::GatherScatterCommunicationTest" );

GatherScatterCommunicationTest:: GatherScatterCommunicationTest ()
:
  TestCase ( "m2n::tests::GatherScatterCommunicationTest" )
{}

void GatherScatterCommunicationTest:: run ()
{
  preciceTrace ( "run" );
# ifndef PRECICE_NO_MPI
  typedef utils::Parallel Par;
  if (Par::getCommunicatorSize() > 3){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1, 2 , 3;
    MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
    if (Par::getProcessRank() <= 3){
      Par::setGlobalCommunicator(comm);
      testMethod ( testSendReceiveAll );
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
# endif // not PRECICE_NO_MPI
}

void GatherScatterCommunicationTest:: testSendReceiveAll ()
{
  preciceTrace ( "testSendReceiveAll" );
  assertion ( utils::Parallel::getCommunicatorSize() == 4 );

  com::PtrCommunication participantCom = com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::PtrGlobalCommunication globalCom = m2n::PtrGlobalCommunication(new m2n::GatherScatterCommunication(participantCom));
  com::PtrCommunication masterSlaveCom = com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0){ //Participant 1
    utils::Parallel::initialize ( NULL, NULL, "Part1" );
    globalCom->acceptConnection ( "Part1", "Part2Master", 0, 1);
  }
  else if(utils::Parallel::getProcessRank() == 1){//Participant 2 - Master
    utils::Parallel::initialize ( NULL, NULL, "Part2Master" );
    globalCom->requestConnection ( "Part1", "Part2Master", 0, 1);
  }
  else if(utils::Parallel::getProcessRank() == 2){//Participant 2 - Slave1
    utils::Parallel::initialize ( NULL, NULL, "Part2Slaves");
  }
  else if(utils::Parallel::getProcessRank() == 3){//Participant 2 - Slave2
    utils::Parallel::initialize ( NULL, NULL, "Part2Slaves");
  }

  if(utils::Parallel::getProcessRank() == 1){//Master
    masterSlaveCom->acceptConnection ( "Part2Master", "Part2Slaves", 0, 1);
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    masterSlaveCom->requestConnection( "Part2Master", "Part2Slaves", 0, 2 );
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    masterSlaveCom->requestConnection( "Part2Master", "Part2Slaves", 1, 2 );
  }


  int dimensions = 2;
  int numberOfVertices = 6;
  bool flipNormals = false;
  int valueDimension = 1;
  utils::DynVector offset ( dimensions, 0.0 );

  if (utils::Parallel::getProcessRank() == 0){ //Part1
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = false;
    mesh::PtrMesh pMesh(new mesh::Mesh("Mesh", dimensions, flipNormals));
    utils::DynVector values(numberOfVertices);
    assignList(values) = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
    globalCom->sendAll(&values,numberOfVertices,0,pMesh,valueDimension);
    globalCom->receiveAll(&values,numberOfVertices,0,pMesh,valueDimension);
    validate(values[0]==2.0);
    validate(values[1]==4.0);
    validate(values[2]==6.0);
    validate(values[3]==16.0);
    validate(values[4]==10.0);
    validate(values[5]==12.0);

  }
  else{
    mesh::PtrMesh pMesh(new mesh::Mesh("Mesh", dimensions, flipNormals));

    if(utils::Parallel::getProcessRank() == 1){//Master
      utils::MasterSlave::_rank = 0;
      utils::MasterSlave::_size = 3;
      utils::MasterSlave::_slaveMode = false;
      utils::MasterSlave::_masterMode = true;

      pMesh->setGlobalNumberOfVertices(numberOfVertices);
      pMesh->getVertexDistribution()[0].push_back(0);
      pMesh->getVertexDistribution()[0].push_back(1);
      pMesh->getVertexDistribution()[0].push_back(3);
      pMesh->getVertexDistribution()[2].push_back(2);
      pMesh->getVertexDistribution()[2].push_back(3);
      pMesh->getVertexDistribution()[2].push_back(4);
      pMesh->getVertexDistribution()[2].push_back(5);

      utils::DynVector values(3);
      assignList(values) = 0.0, 0.0, 0.0;
      globalCom->receiveAll(&values,3,0,pMesh,valueDimension);
      validate(values[0]==1.0);
      validate(values[1]==2.0);
      validate(values[2]==4.0);
      values = values * 2;
      globalCom->sendAll(&values,3,0,pMesh,valueDimension);
    }
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      utils::MasterSlave::_rank = 1;
      utils::MasterSlave::_size = 3;
      utils::MasterSlave::_slaveMode = true;
      utils::MasterSlave::_masterMode = false;
      utils::DynVector values(0);
      globalCom->receiveAll(&values,0,0,pMesh,valueDimension);
      globalCom->sendAll(&values,0,0,pMesh,valueDimension);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      utils::MasterSlave::_rank = 2;
      utils::MasterSlave::_size = 3;
      utils::MasterSlave::_slaveMode = true;
      utils::MasterSlave::_masterMode = false;
      utils::DynVector values(4);
      assignList(values) = 0.0, 0.0, 0.0, 0.0;
      globalCom->receiveAll(&values,4,0,pMesh,valueDimension);
      validate(values[0]==3.0);
      validate(values[1]==4.0);
      validate(values[2]==5.0);
      validate(values[3]==6.0);
      values = values * 2;
      globalCom->sendAll(&values,4,0,pMesh,valueDimension);
    }

  }
  utils::MasterSlave::_slaveMode = false;
  utils::MasterSlave::_masterMode = false;
  utils::Parallel::synchronizeProcesses();
}

}}} // namespace precice, m2n, tests

#endif // PRECICE_NO_MPI
