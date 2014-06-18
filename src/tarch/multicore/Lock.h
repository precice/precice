// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_MULTICORE_LOCK_H_
#define _TARCH_MULTICORE_LOCK_H_


namespace tarch {
  namespace multicore {
    class BooleanSemaphore;
    class Lock;
  }
}


class tarch::multicore::Lock {
  private:
    BooleanSemaphore&  _semaphore;
    bool               _lockIsAquired;
  public:
    Lock( tarch::multicore::BooleanSemaphore& semaphore, bool aquireLockImmediately = true );
    ~Lock();

    void lock();
    void free();
};

#endif
