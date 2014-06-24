#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/logging/Log.h"


#include <limits>


#include <tbb/tbb_machine.h>
#include <tbb/task.h>


int        tarch::multicore::BooleanSemaphore::_pauseCounter(1);
const int  tarch::multicore::BooleanSemaphore::_pauseBeforeYield(32);
const int  tarch::multicore::BooleanSemaphore::_counterThresholdForWarning(std::numeric_limits<int>::max() - 20);


tarch::multicore::BooleanSemaphore::BooleanSemaphore() {
}


tarch::multicore::BooleanSemaphore::~BooleanSemaphore() {
}


void tarch::multicore::BooleanSemaphore::enterCriticalSection() {
  // This creates a seg fault
//  bool gotLock = _lock.try_lock();
//  if (!gotLock) {
//    while (!gotLock) {
//      sendCurrentTaskToBack("enterCriticalSection()");
//      gotLock = _lock.try_lock();
//    }
//    continueWithTask();
//  }
  _lock.lock();
}


void tarch::multicore::BooleanSemaphore::leaveCriticalSection() {
  _lock.unlock();
}


void tarch::multicore::BooleanSemaphore::sendCurrentTaskToBack(const std::string& methodTrace) {
  static tarch::logging::Log  _log( "tarch::multicore::BooleanSemaphore" );
  if (_pauseCounter < _pauseBeforeYield) {
    __TBB_Pause(_pauseCounter);
    _pauseCounter*=2;
  }
  else {
    if (_pauseCounter>_counterThresholdForWarning && _pauseCounter != std::numeric_limits<int>::max()) {
      _pauseCounter = std::numeric_limits<int>::max();
      logWarning( "sendCurrentTaskToBack(string)", "probably running into deadlock or inefficient behaviour in " << methodTrace );
    }
    else {
      _pauseCounter++;
    }
    __TBB_Yield();
  }
}


void tarch::multicore::BooleanSemaphore::continueWithTask() {
  _pauseCounter = 1;
}
