// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MPIDirectCommunicationTest.hpp"
#include "../MPIDirectCommunication.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::com::tests::MPIDirectCommunicationTest)

namespace precice {
namespace com {
namespace tests {

logging::Logger MPIDirectCommunicationTest::
  _log ("precice::com::tests::MPIDirectCommunicationTest");

MPIDirectCommunicationTest:: MPIDirectCommunicationTest ()
:
  TestCase ("precice::com::tests::MPIDirectCommunicationTest")
{}

void MPIDirectCommunicationTest:: run ()
{
# ifndef PRECICE_NO_MPI
  if ( utils::Parallel::getCommunicatorSize() > 1 ) {
    testMethod ( testSendReceiveTwoProcesses );
  }
  if ( utils::Parallel::getCommunicatorSize() > 2 ) {
    testMethod ( testSendReceiveThreeProcesses );
  }
# endif
}

#ifndef PRECICE_NO_MPI

void MPIDirectCommunicationTest:: testSendReceiveTwoProcesses ()
{
  preciceTrace ( "testSendReceiveTwoProcesses()" );
  typedef utils::Parallel Par;
  Par::synchronizeProcesses();

  std::vector<int> ranks;
  ranks += 0, 1;
  MPI_Comm comm = Par::getRestrictedCommunicator ( ranks );
  if ( Par::getProcessRank() < 2 ) {
    Par::setGlobalCommunicator(comm);
    validateEquals ( Par::getCommunicatorSize(), 2 );
    MPIDirectCommunication communication;
    std::string nameEven ( "even" );
    std::string nameOdd  ( "odd" );

    if ( Par::getProcessRank() == 0 ) {
      Par::splitCommunicator( nameEven );
      communication.acceptConnection ( nameEven, nameOdd, 0, 1 );
      int message = 1;
      communication.send ( message, 0 );
      communication.receive ( message, 0 );
      validateEquals ( message, 2 );
      communication.closeConnection();
    }
    else if ( Par::getProcessRank() == 1 ) {
      Par::splitCommunicator( nameOdd );
      communication.requestConnection ( nameEven, nameOdd, 0, 1 );
      int message = -1;
      communication.receive ( message, 0 );
      validateEquals ( message, 1 );
      message = 2;
      communication.send ( message, 0 );
      communication.closeConnection();
    }
    Par::setGlobalCommunicator(Par::getCommunicatorWorld());
  }
}

void MPIDirectCommunicationTest:: testSendReceiveThreeProcesses ()
{
  preciceTrace ( "testSendReceiveThreeProcesses()" );
  utils::Parallel::synchronizeProcesses ();

  std::vector<int> ranks;
  ranks += 0, 1, 2;
  MPI_Comm comm = utils::Parallel::getRestrictedCommunicator ( ranks );
  if ( utils::Parallel::getProcessRank() < 3 ) {
    utils::Parallel::setGlobalCommunicator(comm);
    std::string nameProcess0 ( "one" );
    std::string nameProcess1 ( "two" );
    std::string nameProcess2 ( "three" );

    utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
  }
}

#endif // not PRECICE_NO_MPI



}}} // namespace precice, com, tests
