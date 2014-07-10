// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_MULTICORE_BOOLEAN_SEMAPHORE_TBB_H_
#define _TARCH_MULTICORE_BOOLEAN_SEMAPHORE_TBB_H_


namespace tarch {
  namespace multicore {
    class BooleanSemaphore;
    class Lock;
  }
}



#if defined(SharedTBB)

#include <string>
#include <tbb/spin_mutex.h>


class tarch::multicore::BooleanSemaphore {
  private:
    friend class tarch::multicore::Lock;

    static int                            _pauseCounter;
    static const int                      _pauseBeforeYield;
    static const int                      _counterThresholdForWarning;

    tbb::spin_mutex _lock;

    /**
     * Waits until I can enter the critical section.
     */
    void enterCriticalSection();

    /**
     * Tells the semaphore that it is about to leave.
     */
    void leaveCriticalSection();

    /**
     * You may not copy a semaphore
     */
    BooleanSemaphore( const BooleanSemaphore& semaphore ) {}

    /**
     * You may not copy a semaphore
     */
    BooleanSemaphore& operator=( const BooleanSemaphore& semaphore ) {return *this;}
  public:
    BooleanSemaphore();
    ~BooleanSemaphore();

    /**
     * Send task to background
     *
     * The operation first analyses the _pauseCounter reset by
     * continueWithTask().
     *
     * - If it is smaller than _pauseBeforeYield, we double it, pause, and
     *   return. The pause is realised due to __TBB_Pause which implies a
     *   couple of nops. A thread hence is not suspended, it keeps on running
     *   but it does not disturb other processes. If only this variant would
     *   exist, tasks could starve.
     * - If however the _pauseCounter exceeds, the operation invokes
     *   __TBB_Yield. Furthermore, it writes an error message once.
     */
    static void sendCurrentTaskToBack(const std::string& methodTrace);

    /**
     * Each sendCurrentTaskToBack() should be followed by a continueWithTask().
     */
    static void continueWithTask();

};
#endif


#endif
