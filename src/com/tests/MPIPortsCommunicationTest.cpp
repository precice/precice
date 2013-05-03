// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MPIPortsCommunicationTest.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::com::tests::MPIPortsCommunicationTest)

namespace precice {
namespace com {
namespace tests {

tarch::logging::Log MPIPortsCommunicationTest::
  _log ("precice::com::tests::MPIPortsCommunicationTest");

MPIPortsCommunicationTest:: MPIPortsCommunicationTest()
:
  TestCase ("precice::com::tests::MPIPortsCommunicationTest")
{}

void MPIPortsCommunicationTest:: run()
{
# ifndef PRECICE_NO_MPI
  if ( utils::Parallel::getCommunicatorSize() > 1 ){
    testMethod ( testSendReceiveTwoProcesses );
  }
# endif
}

#ifndef PRECICE_NO_MPI

void MPIPortsCommunicationTest:: testSendReceiveTwoProcesses()
{
  preciceTrace ( "testSendReceiveTwoProcesses()" );
  typedef utils::Parallel Par;
  Par::synchronizeProcesses();

  std::vector<int> ranks;
  ranks += 0, 1;
  MPI_Comm comm = Par::getRestrictedCommunicator(ranks);
  if ( Par::getProcessRank() < 2 ) {
    Par::setGlobalCommunicator(comm);
    validateEquals ( Par::getCommunicatorSize(), 2 );
    MPIPortsCommunication communication("");
    std::string nameEven ( "even" );
    std::string nameOdd  ( "odd" );

    if ( Par::getProcessRank() == 0 ) {
      //Par::initialize ( NULL, NULL, nameEven );
      communication.acceptConnection ( nameEven, nameOdd, 0, 1 );
      int message = 1;
      communication.send ( message, 0 );
      communication.receive ( message, 0 );
      validateEquals ( message, 2 );
      communication.closeConnection ();
    }
    else if ( Par::getProcessRank() == 1 ) {
      //Par::initialize ( NULL, NULL, nameOdd );
      communication.requestConnection ( nameEven, nameOdd, 0, 1 );
      int message = -1;
      communication.receive ( message, 0 );
      validateEquals ( message, 1 );
      message = 2;
      communication.send ( message, 0 );
      communication.closeConnection ();
    }
    Par::setGlobalCommunicator(Par::getCommunicatorWorld());
  }
}

#endif // not PRECICE_NO_MPI

}}} // namespace precice, com, tests
