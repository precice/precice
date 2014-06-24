#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/multicore/openMP/Core.h"
#include "tarch/Assertions.h"

#include <iostream>

#ifdef SharedOMP
#include <omp.h>
#endif

tarch::multicore::openMP::Core::Core():
  _numberOfThreads(omp_get_num_procs()) {
  omp_set_num_threads(_numberOfThreads);
}


tarch::multicore::openMP::Core::~Core() {

}


tarch::multicore::openMP::Core& tarch::multicore::openMP::Core::getInstance() {
  static tarch::multicore::openMP::Core singleton;
  return singleton;
}


void tarch::multicore::openMP::Core::shutDown() {
  _numberOfThreads = -1;
}


void tarch::multicore::openMP::Core::configure( int numberOfThreads ) {
  if (numberOfThreads==0) {
    _numberOfThreads = omp_get_num_procs();
  }
  else {
    _numberOfThreads = numberOfThreads;
    omp_set_num_threads(_numberOfThreads);
  }
}


int tarch::multicore::openMP::Core::getNumberOfThreads() const {
  return _numberOfThreads;
}

bool tarch::multicore::openMP::Core::isInitialised() const {
  return _numberOfThreads > 0;
}
