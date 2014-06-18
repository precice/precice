#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/logging/Log.h"
#include "tarch/logging/CommandLineLogger.h"
#include "tarch/logging/CCALogger.h"
#include "tarch/compiler/CompilerSpecificSettings.h"
#include "tarch/Assertions.h"

/**
 * For the machine name. If it doesn't work, switch it off in the file
 * CompilerSpecificSettings.h.
 */
#ifdef CompilerHasUTSName
#include <sys/utsname.h>
#endif

#ifdef Parallel
#include "tarch/parallel/Node.h"
#endif

#include <time.h>


tarch::logging::Log::Log(const std::string& className):
_className( className ) {
}


tarch::logging::Log::~Log() {
}


#if defined(Debug) && !defined(LogOff)
void tarch::logging::Log::debug(const std::string& methodName, const std::string& message) const {
  UsedLogService::getInstance().debug(getTimeStampMS(),getTimeStampHumanReadable(),getMachineInformation(),getTraceInformation(methodName),message);
}

void tarch::logging::Log::debugMasterOnly(const std::string& methodName, const std::string& message) const {
#ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
#endif
    debug(methodName, message);
#ifdef Parallel
  }
#endif
}

#endif


#if !defined(LogOff)
void tarch::logging::Log::info(const std::string& methodName, const std::string& message) const {
  UsedLogService::getInstance().info(getTimeStampMS(),getTimeStampHumanReadable(),getMachineInformation(),getTraceInformation(methodName),message);
}

void tarch::logging::Log::infoMasterOnly(const std::string& methodName, const std::string& message) const {
#ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
#endif
    info(methodName, message);
#ifdef Parallel
  }
#endif
}

void tarch::logging::Log::warning(const std::string& methodName, const std::string& message) const {
  UsedLogService::getInstance().warning(getTimeStampMS(),getTimeStampHumanReadable(),getMachineInformation(),getTraceInformation(methodName),message);
}

void tarch::logging::Log::warningMasterOnly(const std::string& methodName, const std::string& message) const {
#ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
#endif
    warning(methodName, message);
#ifdef Parallel
  }
#endif
}

void tarch::logging::Log::error(const std::string& methodName, const std::string& message) const {
  UsedLogService::getInstance().error(getTimeStampMS(),getTimeStampHumanReadable(),getMachineInformation(),getTraceInformation(methodName),message);
}

void tarch::logging::Log::indent( bool indent, const std::string& trace, const std::string& message ) const {
#if defined(Debug)
  UsedLogService::getInstance().indent( indent, trace, message );
#endif
}
#endif


std::string tarch::logging::Log::getMachineInformation() const {
  std::ostringstream message;

#ifdef CompilerHasUTSName
  utsname* utsdata = new utsname();
  assertion( utsdata!=NULL );
  uname(utsdata);

  message << "[" << utsdata->nodename << "]";

#ifdef Parallel
  message << ",";
#endif

  delete utsdata;
#endif

#ifdef Parallel
  if (tarch::parallel::Node::getInstance().isInitialised()) {
    message << "rank:" << tarch::parallel::Node::getInstance().getRank();
  }
  else {
    message << "rank:not-initialised-yet";
  }
#endif

  return message.str();
}


long int tarch::logging::Log::getTimeStampMS() const {
  time_t* timeStamp = new time_t();
  assertion( timeStamp!=NULL );
  time(timeStamp);
  long int result = *timeStamp;
  delete timeStamp;
  return result;
}


std::string tarch::logging::Log::getTimeStampHumanReadable() const {
  // calender time: create struct and get time from system
  time_t* timeStamp = new time_t();
  assertion( timeStamp!=NULL );
  time(timeStamp);

  // Break down time into hour, seconds, ...
  // Note that time is only a substructure of timeStamp. Therefore the pointer
  // to time may not be deleted.
  tm*     time      = localtime(timeStamp);
  assertion( time!=NULL );

  std::ostringstream message;

  // write all information
  if (time->tm_hour<10) {
    message << "0";
  }
  message << time->tm_hour << ":";

  if (time->tm_min<10) {
    message << "0";
  }
  message << time->tm_min << ":";

  if (time->tm_sec<10) {
    message << "0";
  }
  message << time->tm_sec;

  delete timeStamp;

  return message.str();
}


std::string tarch::logging::Log::getTraceInformation( const std::string& methodName ) const {
  return _className + "::" + methodName;
}
