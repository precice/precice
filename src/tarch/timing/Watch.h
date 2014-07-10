// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_TIMING_WATCH_H_
#define _TARCH_TIMING_WATCH_H_


#include <ctime>
#include "tarch/logging/Log.h"

#ifdef __APPLE__
#include <mach/mach_time.h>
#include <mach/mach.h>
#include <mach/clock.h>
#endif

namespace tarch {
  namespace timing {
    class Watch;
  }
}



#ifdef SharedOMP
#include <omp.h>
#elif SharedTBB
#include <tbb/tick_count.h>
#endif


/**
 * A simple class that has to be included to measure the clock ticks required
 * for an operation. To use it you've to create an instance of this class at
 * the beginning of the operation block you want to measure the time spent
 * within it:
 *
 * <pre>
 *   void anyOperation() {
 *     utils::Watch watch("MyClass","anyOperation()");
 *     ...
 *   }
 * </pre>
 *
 * The result the of the measurement is written to log.info level. Note that
 * this operation works within operation blocks (e.g. a for-loop), too.
 *
 * For an extended measurement of calender (i.e. user) time and a saving of the
 * data and access via getters, we implemented parts of what Markus Langlotz
 * gave us for the ancient stokes project. Additionally, the counter overflow
 * and the non-access to the measured time were reasons for introducing this
 * extended version. Furthermore, we now may choose specific computation
 * intervals to be counted in between the start and the stop command.
 *
 * @version $Revision: 1.9 $
 * @author  Tobias Weinzierl, Tobias Neckel
 */
class tarch::timing::Watch {
  private:
    /**
     * Log device the result is written to.
     */
    tarch::logging::Log     _log;

    /**
     * Flag to distinguish the (original) standard watch. Is used in the
     * destructor.
     */
    bool _plotResultInDestructor;

    /**
     * Flag to check if only statistics should be output in the parallel case.
     * only relevant in the parallel case
     */
    bool _onlyOutputAverages;

    /**
     * Stores the name of the operation the watch is used within.
     */
    std::string    _operationName;

    /**
     * Holds the clock ticks at the beginning of the time measurement.
     */
    std::clock_t   _startClockTicks;

    #ifdef SharedOMP
    double          _startTime;
    #elif SharedTBB
    tbb::tick_count _startTime;
    #else
    /**
     * Holds the time at the beginning of the time measurement.
     */
    double          _startTime;
    //std::time_t    _startTime;
    #endif

    /**
     * Holds the elapsed processor time.
     */
    std::clock_t   _elapsedClockTicks;

    /**
     * Holds the elapsed calendar time.
     */
    double         _elapsedTime;

    /**
     * Has stopTimer() been called before.
     */
    bool _calledStopManually;

    #ifdef __APPLE__
    clock_serv_t cclock;
    #endif

  public:
    /**
     * Construct a watch
     *
     * @param className     Name of the class the watch is used within (for log).
     * @param operationName Name of the operation the watch is used within.
     * @param extendedWatchType Bool for extended.
     */
    Watch(
      const std::string& className,
      const std::string& operationName,
      const bool         plotResultInDestructor,
      const bool         onlyOutputAverages = false
    );

    /**
     * For standard version (): Stops the watch and plots the time spent. For
     * the extended version, this is the usual destructor without any special
     * actions.
     *
     * In the parallel one can choose to only some general statistics and not
     * all timings of all processes, namely, minimum, maximum and average value.
     * This functionality is activated by setting _onlyOutputAverages true.
     */
    virtual ~Watch();

    /**
     * (Re)Start the Timer
     *
     * This method starts the timer. Actually, it restarts the time as the
     * constructor also invokes this operation.
     */
    void startTimer();


    /**
     * This method stops the timer
     */
    void stopTimer();

    /**
     * Return CPU Time in Seconds
     *
     * This method returns the elapsed cpu time between the start and stop
     * command of the timer, i.e. the clock ticks actually spent by the process.
     * Take care: This counter might overflow, especially on a 32 bit
     * architecture.
     *
     * !!! Multithreading
     *
     * If you use multithreading, cpu time sums up the clock ticks invested by
     * the individual cores, i.e. if you switch from a one core to a dual core
     * machine, the cpu time should roughly double.
     *
     * @return CPU time in seconds
     */
    double getCPUTime();

    /**
     * Equals getCPUTime() but returns the clock ticks instead of the time in
     * seconds.
     */
    std::clock_t getCPUTicks();

    /**
     * This method returns the elapsed calendar time between the start and stop
     * command of the timer, i.e. the real world time. The result is given in
     * seconds.
     */
    double getCalendarTime();
};

#endif
