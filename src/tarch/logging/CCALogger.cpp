#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#if defined(CCA)
#include "tarch/logging/CCALogger.h"
#include "tarch/logging/cca/Log.h"
#include "tarch/logging/CommandLineLogger.h"


#include "peano/integration/cca/UsesPortService.h"


#include "tarch/multicore/Lock.h"

#include "tarch/Assertions.h"

#ifdef Parallel
#include "tarch/parallel/Node.h"
#endif

#include <iostream>


const std::string tarch::logging::CCALogger::CCAUsesPortLog( "tarch::logging::cca::Log" );


tarch::logging::CCALogger::CCALogger():
  _isCalled(false) {
}


tarch::logging::CCALogger::~CCALogger() {
  assertion( !_isCalled );
}


tarch::logging::CCALogger& tarch::logging::CCALogger::getInstance() {
  static tarch::logging::CCALogger singleton;
  return singleton;
}


void tarch::logging::CCALogger::debug(const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message) {
  if (!_isCalled) {
    _isCalled = true;

    tarch::multicore::Lock lock( _semaphore );

    #ifdef Parallel
    const int rank = tarch::parallel::Node::getInstance().getRank();
    #else
    const int rank = -1;
    #endif

    if (peano::integration::cca::UsesPortService::getInstance().hasUsesPort( CCAUsesPortLog )) {
      peano::integration::cca::UsesPortService::getInstance().getUsesPort<tarch::logging::cca::Log>(CCAUsesPortLog).debug(timestampMS,timestampHumanReadable,rank,machineName,trace,message);
    }
    else {
      CommandLineLogger::getInstance().debug(timestampMS,timestampHumanReadable,machineName,trace,message);
    }

    _isCalled = false;
  }
  else {
    CommandLineLogger::getInstance().error(timestampMS,timestampHumanReadable,machineName,trace,"recursive try to write debug statement: " + message);
  }
}


void tarch::logging::CCALogger::info(const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message) {
  if (!_isCalled) {
    _isCalled = true;

    tarch::multicore::Lock lock( _semaphore );

    #ifdef Parallel
    const int rank = tarch::parallel::Node::getInstance().getRank();
    #else
    const int rank = -1;
    #endif

    if (peano::integration::cca::UsesPortService::getInstance().hasUsesPort( CCAUsesPortLog )) {
      peano::integration::cca::UsesPortService::getInstance().getUsesPort<tarch::logging::cca::Log>(CCAUsesPortLog).info(timestampMS,timestampHumanReadable,rank,machineName,trace,message);
    }
    else {
      CommandLineLogger::getInstance().info(timestampMS,timestampHumanReadable,machineName,trace,message);
    }

    _isCalled = false;
  }
  else {
    CommandLineLogger::getInstance().error(timestampMS,timestampHumanReadable,machineName,trace,"recursive try to write info statement: " + message);
  }
}


void tarch::logging::CCALogger::warning(const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message) {
  if (!_isCalled) {
   _isCalled = true;

    tarch::multicore::Lock lock( _semaphore );

    #ifdef Parallel
    const int rank = tarch::parallel::Node::getInstance().getRank();
    #else
    const int rank = -1;
    #endif

    if (peano::integration::cca::UsesPortService::getInstance().hasUsesPort( CCAUsesPortLog )) {
      peano::integration::cca::UsesPortService::getInstance().getUsesPort<tarch::logging::cca::Log>(CCAUsesPortLog).warning(timestampMS,timestampHumanReadable,rank,machineName,trace,message);
    }

    CommandLineLogger::getInstance().warning(timestampMS,timestampHumanReadable,machineName,trace,message);

    _isCalled = false;
  }
  else {
    CommandLineLogger::getInstance().error(timestampMS,timestampHumanReadable,machineName,trace,"recursive try to write warning statement: " + message);
  }
}


void tarch::logging::CCALogger::error(const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message) {
  if (!_isCalled) {
    _isCalled = true;

    tarch::multicore::Lock lock( _semaphore );

    #ifdef Parallel
    const int rank = tarch::parallel::Node::getInstance().getRank();
    #else
    const int rank = -1;
    #endif

    if (peano::integration::cca::UsesPortService::getInstance().hasUsesPort( CCAUsesPortLog )) {
      peano::integration::cca::UsesPortService::getInstance().getUsesPort<tarch::logging::cca::Log>(CCAUsesPortLog).error(timestampMS,timestampHumanReadable,rank,machineName,trace,message);
    }

    CommandLineLogger::getInstance().error(timestampMS,timestampHumanReadable,machineName,trace,message);

    _isCalled = false;
  }
  else {
    CommandLineLogger::getInstance().error(timestampMS,timestampHumanReadable,machineName,trace,"recursive try to write error statement: " + message);
  }
}


void tarch::logging::CCALogger::indent( bool indent ) {
  CommandLineLogger::getInstance().indent( indent );
}

#endif
