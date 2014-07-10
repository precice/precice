#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/multicore/openMP/BooleanSemaphore.h"

#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/multicore/openMP/Core.h"

tarch::multicore::BooleanSemaphore::BooleanSemaphore() {
  omp_init_lock(&_lock);
}


tarch::multicore::BooleanSemaphore::~BooleanSemaphore() {
  #ifdef CompilerICC
  // Fix for Intel Compiler due to a bug in the OpenMP support of the icc (see http://software.intel.com/en-us/forums/showthread.php?t=72204)
  if(tarch::multicore::openMP::Core::getInstance().isInitialised()) {
    omp_destroy_lock(&_lock);
  }
  #else
  omp_destroy_lock(&_lock);
  #endif
}


void tarch::multicore::BooleanSemaphore::enterCriticalSection() {
  omp_set_lock(&_lock);
}


void tarch::multicore::BooleanSemaphore::leaveCriticalSection() {
  omp_unset_lock(&_lock);
}



void tarch::multicore::BooleanSemaphore::sendCurrentTaskToBack(const std::string& methodTrace) {
}


void tarch::multicore::BooleanSemaphore::continueWithTask() {
}
