#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/timing/Watch.h"
#include <time.h>

#include <sstream>

#ifdef Parallel
#include "tarch/parallel/Statistics.h"
#endif

#if !defined(SharedTBB) && !defined(SharedOMP)
/**
 * Forward declaration
 *
 * We need that dummy if and only if we are not doing shared memory stuff. If
 * we do shared memory, we rely on the shared memory libs' operations to get
 * the right time instead of this timespec struct.
 */
struct timespec ts;
#endif


tarch::timing::Watch::Watch(
  const std::string& className,
  const std::string& operationName,
  const bool         plotResultInDestructor,
  const bool         onlyOutputAverages
  ):
  _log( className ),
  _plotResultInDestructor(plotResultInDestructor),
  _onlyOutputAverages(onlyOutputAverages),
  _operationName( operationName ),
  _startClockTicks( 0 ),
  _startTime(),
  _elapsedClockTicks( 0 ),
  _elapsedTime( 0 ),
  _calledStopManually(false)
  {
  #ifdef __APPLE__
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  #endif
  startTimer();
}


tarch::timing::Watch::~Watch() {
  if (!_calledStopManually) {
    stopTimer();
  }

  if (_plotResultInDestructor) {
    std::ostringstream message;
  #ifdef Parallel
    if(!_onlyOutputAverages) {
  #endif
      // standard serial output
      message << "total number of clock ticks within block (cpu-time,calendar-time): "
          << "(" << getCPUTime() << "s"
          << "," << getCalendarTime() << "s"
          << ")";
      _log.info( _operationName, message.str() );
  #ifdef Parallel
    } else {
      if(tarch::parallel::Node::getInstance().isGlobalMaster())
        message << "process statistics of clock ticks within block \ncpu-time:";
      tarch::parallel::Statistics::addStatisticsToStream(getCPUTime(),"s",message);
      if(tarch::parallel::Node::getInstance().isGlobalMaster())
        message <<"\ncalendar-time: ";
      tarch::parallel::Statistics::addStatisticsToStream(getCalendarTime(),"s",message);
      _log.infoMasterOnly( _operationName, message.str() );
    }
  #endif
  }

  #ifdef __APPLE__
  mach_port_deallocate(mach_task_self(), cclock);
  #endif
}


void tarch::timing::Watch::startTimer() {
  _startClockTicks = std::clock();
  #ifdef SharedOMP
  _startTime       = omp_get_wtime();
  #elif SharedTBB
  _startTime       = tbb::tick_count::now();
  #elif defined(__APPLE__)
  mach_timespec_t mts;
  clock_get_time(cclock, &mts);
  _startTime = (double)mts.tv_sec + (double)mts.tv_nsec * 1e-09;
  #else
   if( clock_gettime(CLOCK_REALTIME, &ts) == 0 ) {
	   _startTime = (double)ts.tv_sec + (double)ts.tv_nsec * 1e-09;
   }
  #endif
}


void tarch::timing::Watch::stopTimer() {
  _elapsedClockTicks   = clock() - _startClockTicks;
  #ifdef SharedOMP
  _elapsedTime         = omp_get_wtime() - _startTime;
  #elif SharedTBB
  _elapsedTime         = (tbb::tick_count::now() - _startTime).seconds();
  #elif defined(__APPLE__)
  mach_timespec_t mts;
  clock_get_time(cclock, &mts);
  _startTime = (double)mts.tv_sec + (double)mts.tv_nsec * 1e-09;
  #else
  clock_gettime(CLOCK_REALTIME, &ts);
  _elapsedTime         = ((double)ts.tv_sec + (double)ts.tv_nsec * 1e-09)-_startTime;
  //difftime ( std::time(NULL), _startTime );
  #endif
  _calledStopManually  = true;
}


double tarch::timing::Watch::getCPUTime() {
  double lhs = static_cast<double>(_elapsedClockTicks);
  double rhs = static_cast<double>(CLOCKS_PER_SEC);
  return lhs/rhs;
}


std::clock_t tarch::timing::Watch::getCPUTicks() {
  return _elapsedClockTicks;
}


double tarch::timing::Watch::getCalendarTime()  {
  return _elapsedTime;
}
