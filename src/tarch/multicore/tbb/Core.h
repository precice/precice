// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_MULTICORE_TBB_CORE_H_
#define _TARCH_MULTICORE_TBB_CORE_H_


#include <tbb/task_scheduler_init.h>

namespace tarch {
  namespace multicore {
    namespace tbb {
      class Core;
    }
  }
}


/**
 * TBB Core
 *
 * Abstraction of the TBB routines. This class is a singleton.
 *
 * @author Tobias Weinzierl
 */
class tarch::multicore::tbb::Core {
  private:
    Core();

    int                         _numberOfThreads;
    ::tbb::task_scheduler_init  _task_scheduler_init;
  public:
    /**
     * Destructor
     */
    ~Core();

    /**
     * @return Singleton instance
     */
    static Core& getInstance();

    /**
     * Configure the whole thing. If numberOfThreads equals 0, the core is
     * using the number of standard threads.
     *
     *
     * @param numberOfThreads Numer of threads that shall be used. This
     *        parameter either is greater than zero (which defines the number
     *        of threads) or it equals zero which means that the code should
     *        use the default number of threads. The latter value fits to the
     *        semantics of
     *        tarch::multicore::configurations::CoreConfiguration::getNumberOfThreads()
     */
    void configure( int numberOfThreads );

    /**
     * Shutdown parallel environment.
     */
    void shutDown();

    /**
     * @return Shared memory environment is up and runnning.
     */
    bool isInitialised() const;

    /**
     * Returns the number of threads that is used by TBB. This routine usually
     * is not of interest at all as TBB should do all the thread management
     * automatically. However, the user interface plots some information on the
     * number of threads used, and sometimes I found it useful.
     *
     * @return Number of threads available.
     */
    int getNumberOfThreads() const;
};


#endif
