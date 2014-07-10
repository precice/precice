// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_MULTICORE_BOOLEAN_SEMAPHORE_COBRA_H_
#define _TARCH_MULTICORE_BOOLEAN_SEMAPHORE_COBRA_H_



namespace tarch {
  namespace multicore {
    class BooleanSemaphore;
    class Lock;
  }
}



#if defined(SharedCobra)

#include <string>
#include <cobra/thread.hpp>


class tarch::multicore::BooleanSemaphore {
  private:
    friend class tarch::multicore::Lock;

    static int              _pauseCounter;
    static const int        _pauseBeforeYield;
    static const int        _counterThresholdForWarning;

    ::cobra::mutex          _mutex;

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
     * @todo Currently not supported, i.e. application might starve.
     */
    static void sendCurrentTaskToBack(const std::string& methodTrace);

    /**
     * Each sendCurrentTaskToBack() should be followed by a continueWithTask().
     */
    static void continueWithTask();
};
#endif


#endif
