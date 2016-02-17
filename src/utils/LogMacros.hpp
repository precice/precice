#pragma once 

#include "MasterSlave.hpp"

#include "tarch/Assertions.h"

#include <boost/foreach.hpp>

#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

//#define assertion assertion

#define foreach BOOST_FOREACH

#ifndef PRECICE_NO_MPI
#include "utils/Parallel.hpp"
#define PRECICE_PROCESS_RANK_STREAM \
  "(" << precice::utils::Parallel::getProcessRank() << ") "
#else
#define PRECICE_PROCESS_RANK_STREAM ""
#endif

#define precicePrint(stream) \
   { \
      std::ostringstream conv; \
      conv.setf ( std::ios::showpoint ); \
      conv.setf ( std::ios::fixed ); \
      conv << std::setprecision(16); \
      conv << PRECICE_PROCESS_RANK_STREAM; \
      conv << stream; \
      std::cout << conv.str() << std::endl; \
   }

#include "tarch/logging/Log.h"
#define preciceDebug(stream) \
  logDebug(preciceMethodName, PRECICE_PROCESS_RANK_STREAM << stream)
#define preciceWarning(methodname, stream) \
  logWarning(methodname, PRECICE_PROCESS_RANK_STREAM << stream)
#define preciceInfo(methodname, stream) \
  if(not precice::utils::MasterSlave::_slaveMode)       \
    logInfo(methodname, stream)

/**
 * @brief Standard logging device used in macro preciceDebug.
 */
#define PRECICE_LOGGING_DEVICE _log

#ifdef Debug

#include "Tracer.hpp"

#define preciceTrace(methodname) \
  std::string preciceMethodName(methodname); \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, "");

#define preciceTrace1(methodname, var1) \
   std::string preciceMethodName(methodname); \
   std::ostringstream preciceTraceStream; \
   preciceTraceStream << #var1 << "=" << var1; \
   precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define preciceTrace2(methodname, var1, var2) \
  std::string preciceMethodName(methodname); \
  std::ostringstream preciceTraceStream; \
  preciceTraceStream << #var1 << "=" << var1 << ", " << #var2 << "=" << var2; \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define preciceTrace3(methodname, var1, var2, var3) \
  std::string preciceMethodName(methodname); \
  std::ostringstream preciceTraceStream; \
  preciceTraceStream         << #var1 << "=" << var1  \
                     << ", " << #var2 << "=" << var2  \
                     << ", " << #var3 << "=" << var3; \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define preciceTrace4(methodname, var1, var2, var3, var4) \
  std::string preciceMethodName(methodname); \
  std::ostringstream preciceTraceStream; \
  preciceTraceStream         << #var1 << "=" << var1  \
                     << ", " << #var2 << "=" << var2  \
                     << ", " << #var3 << "=" << var3  \
                     << ", " << #var4 << "=" << var4; \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define preciceTrace5(methodname, var1, var2, var3, var4, var5) \
  std::string preciceMethodName(methodname); \
  std::ostringstream preciceTraceStream; \
  preciceTraceStream         << #var1 << "=" << var1  \
                     << ", " << #var2 << "=" << var2  \
                     << ", " << #var3 << "=" << var3  \
                     << ", " << #var4 << "=" << var4  \
                     << ", " << #var5 << "=" << var5; \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define preciceTrace6(methodname, var1, var2, var3, var4, var5, var6) \
  std::string preciceMethodName(methodname); \
  std::ostringstream preciceTraceStream; \
  preciceTraceStream         << #var1 << "=" << var1  \
                     << ", " << #var2 << "=" << var2  \
                     << ", " << #var3 << "=" << var3  \
                     << ", " << #var4 << "=" << var4  \
                     << ", " << #var5 << "=" << var5  \
                     << ", " << #var6 << "=" << var6; \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#else // Debug

//#define preciceDebug(methodname, message)
#define preciceTrace(methodname)
#define preciceTrace1(methodname, var1)
#define preciceTrace2(methodname, var1, var2)
#define preciceTrace3(methodname, var1, var2, var3)
#define preciceTrace4(methodname, var1, var2, var3, var4)
#define preciceTrace5(methodname, var1, var2, var3, var4, var5)
#define preciceTrace6(methodname, var1, var2, var3, var4, var5, var6)

#endif // ! Debug

/**
 * @brief Needed for macros exiting program execution, to empty logging cache.
 */
//#ifdef _UTILS_LOG_H_
//#define PRECICE_CLOSE_LOGGER tarch::logging::Log::criticalAbort();
//#else
//#define PRECICE_CLOSE_LOGGER
//#endif

/**
 * @brief Wrapper for tarch::logging::Log error output with program exit.
 */
#define preciceError(methodname, message) \
   { \
      std::ostringstream conv; \
      conv << PRECICE_PROCESS_RANK_STREAM << " [PRECICE] ERROR: " << message; \
      PRECICE_LOGGING_DEVICE.error (methodname, conv.str()); \
      std::abort(); \
   }

#define preciceCheck(check, methodname, errormessage) \
   if ( !(check) ) { \
      preciceError(methodname, errormessage); \
   }


