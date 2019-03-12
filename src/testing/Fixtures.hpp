#pragma once

#include "Testing.hpp"
#include "utils/Parallel.hpp"
#include "utils/MasterSlave.hpp"
#include "com/Communication.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/M2N.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/SharedPointer.hpp"

namespace precice
{
namespace testing
{

#ifndef PRECICE_NO_MPI

/// Fixture to create and destroy a master communication
/**
 * Many tests with parallel features need a working master communication. This fixture avoid code duplication for such tests.
 */
struct MasterComFixture {
  MasterComFixture()
  {
    utils::MasterSlave::_communication = com::PtrCommunication(new com::MPIDirectCommunication());
    int size = Par::getCommunicatorSize();

    if (utils::Parallel::getProcessRank() == 0){ //Master
      utils::Parallel::splitCommunicator( "Master" );
      utils::MasterSlave::_rank = 0;
      utils::MasterSlave::_size = size;
      utils::MasterSlave::_slaveMode = false;
      utils::MasterSlave::_masterMode = true;
      utils::MasterSlave::_communication->acceptConnection ( "Master", "Slaves", utils::Parallel::getProcessRank() );
      utils::MasterSlave::_communication->setRankOffset(1);
    }
    else {//Slaves
      assertion(utils::Parallel::getProcessRank() > 0 && utils::Parallel::getProcessRank() < size);
      utils::Parallel::splitCommunicator( "Slaves" );
      utils::MasterSlave::_rank = utils::Parallel::getProcessRank();
      utils::MasterSlave::_size = size;
      utils::MasterSlave::_slaveMode = true;
      utils::MasterSlave::_masterMode = false;
      utils::MasterSlave::_communication->requestConnection( "Master", "Slaves", utils::Parallel::getProcessRank()-1 , size-1 );
    }
  }

  ~MasterComFixture()
  {
    utils::MasterSlave::_communication = nullptr;
    utils::MasterSlave::reset();
    utils::Parallel::clearGroups();
  }
};

/// Fixture to create and destroy an m2n communication between two participants
struct M2NFixture {
  m2n::PtrM2N m2n;

  M2NFixture()
  {
    auto participantCom = com::PtrCommunication(new com::MPIDirectCommunication());
    auto distrFactory = m2n::DistributedComFactory::SharedPointer(new m2n::GatherScatterComFactory(participantCom));
    m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

    if (utils::Parallel::getProcessRank() == 0){
      utils::Parallel::splitCommunicator( "ParticipantOne" );
      m2n->acceptMasterConnection ( "ParticipantOne", "ParticipantTwo");
    }
    else if(utils::Parallel::getProcessRank() == 1){//Master
      utils::Parallel::splitCommunicator( "ParticipantTwo" );
      m2n->requestMasterConnection ( "ParticipantOne", "ParticipantTwo" );
    }
  }

  ~M2NFixture()
  {
    utils::Parallel::clearGroups();
  }
};

#endif // PRECICE_NO_MPI

/// Fixture to split two participants such that both can interact in an integration test
struct SplitParticipantsFixture {
  int participantID;

  SplitParticipantsFixture()
  {
    if(utils::Parallel::getProcessRank()<=1){
      utils::Parallel::splitCommunicator( "ParticipantOne" );
      utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
      utils::Parallel::clearGroups(); //This is important, if the testcase uses MPI communication again
      participantID = 1;
    }
    else {
      assertion(utils::Parallel::getProcessRank() > 1 && utils::Parallel::getProcessRank() < 4);
      utils::Parallel::splitCommunicator( "ParticipantTwo" );
      utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
      utils::Parallel::clearGroups();
      participantID = 2;
    }
  }

  ~SplitParticipantsFixture()
  {
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
  }

};

}}
