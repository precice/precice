// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_MULTICORE_COBRA_CORE_H_
#define _TARCH_MULTICORE_COBRA_CORE_H_


#include <cobra/scheduler.hpp>


namespace tarch {
  namespace multicore {
    namespace cobra {
      class Core;
    }
  }
}


/**
 * Core abstraction for Cobra.
 *
 * This class is a singleton and manages basically the Cobra scheduler.
 * Different to OpenMP and TBB, this singleton always has to be initialised.
 *
 * @author Tobias Weinzierl
 */
class tarch::multicore::cobra::Core {
  private:
    Core();

    /**
     * The cobra scheduler unfortunately does not provide a
     * getNumberOfThreads() operation, so I have to store this
     * information explicitly.
     */
    int                   _numberOfThreads;
    ::cobra::scheduler*   _scheduler;
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
     * Configure the whole thing.
     *
     * @todo Es waere schoen, wenn Cobra default Thread-zahlen akzeptieren wuerde
     *
     * @param numberOfThreads Numer of threads that shall be used. This
     *        parameter has to be greater than zero (which defines the number
     *        of threads)
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

    /**
     * Return scheduler
     */
    ::cobra::scheduler&  getScheduler();
};


#endif
