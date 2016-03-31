#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/multicore/cobra/Core.h"
#include "utils/assertion.hpp"


tarch::multicore::cobra::Core::Core():
  _numberOfThreads(-1),
  _scheduler(0) {
}


tarch::multicore::cobra::Core::~Core() {
  if (_scheduler!=0) {
    shutDown();
  }
}


tarch::multicore::cobra::Core& tarch::multicore::cobra::Core::getInstance() {
  static tarch::multicore::cobra::Core singleton;
  return singleton;
}


void tarch::multicore::cobra::Core::shutDown() {
  delete _scheduler;
  _scheduler       = 0;
  _numberOfThreads = -1;
}


void tarch::multicore::cobra::Core::configure( int numberOfThreads ) {
  assertion( numberOfThreads>0 );

  _scheduler       = new ::cobra::scheduler(numberOfThreads);
  _numberOfThreads = numberOfThreads;
}


int tarch::multicore::cobra::Core::getNumberOfThreads() const {
  assertion( isInitialised() );
  return _numberOfThreads;
}


bool tarch::multicore::cobra::Core::isInitialised() const {
  return _scheduler != 0;
}


::cobra::scheduler&  tarch::multicore::cobra::Core::getScheduler() {
  assertion( isInitialised() );
  return *_scheduler;
}
