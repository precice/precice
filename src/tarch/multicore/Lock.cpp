#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/Assertions.h"
#include "tarch/multicore/Lock.h"
#include "tarch/multicore/BooleanSemaphore.h"


tarch::multicore::Lock::Lock( tarch::multicore::BooleanSemaphore& semaphore, bool aquireLockImmediately ):
  _semaphore(semaphore),
  _lockIsAquired(false) {
  if (aquireLockImmediately) {
    lock();
  }
}


tarch::multicore::Lock::~Lock() {
  if (_lockIsAquired) {
    free();
  }
}


void tarch::multicore::Lock::lock() {
  assertion( !_lockIsAquired );

  _semaphore.enterCriticalSection();

  _lockIsAquired = true;
}


void tarch::multicore::Lock::free() {
  assertion( _lockIsAquired );

  _semaphore.leaveCriticalSection();

  _lockIsAquired = false;
}

