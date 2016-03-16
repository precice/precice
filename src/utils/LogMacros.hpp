#pragma once 

#include "MasterSlave.hpp"

#include <string>
#include <sstream>
#include <iomanip>
#include <cstdlib>

//#define assertion assertion


#ifndef PRECICE_NO_MPI
#include "utils/Parallel.hpp"
#define PRECICE_PROCESS_RANK_STREAM \
  "(" << precice::utils::Parallel::getProcessRank() << ") "
#else
#define PRECICE_PROCESS_RANK_STREAM ""
#endif

#define tprecicePrint(stream) \
   { \
      std::ostringstream conv; \
      conv.setf ( std::ios::showpoint ); \
      conv.setf ( std::ios::fixed ); \
      conv << std::setprecision(16); \
      conv << PRECICE_PROCESS_RANK_STREAM; \
      conv << stream; \
      std::cout << conv.str() << std::endl; \
   }

#include "logging/Logger.hpp"
#define tpreciceDebug(stream) \
  logDebug(preciceMethodName, PRECICE_PROCESS_RANK_STREAM << stream)
#define tpreciceWarning(methodname, stream) \
  logWarning(methodname, PRECICE_PROCESS_RANK_STREAM << stream)
#define tpreciceInfo(methodname, stream) \
  if(not precice::utils::MasterSlave::_slaveMode)       \
    logInfo(methodname, stream)

/**
 * @brief Standard logging device used in macro preciceDebug.
 */
#define PRECICE_LOGGING_DEVICE _log

#ifdef Debug

#include "Tracer.hpp"

#define tpreciceTrace(methodname) \
  std::string preciceMethodName(methodname); \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, "");

#define tpreciceTrace1(methodname, var1) \
   std::string preciceMethodName(methodname); \
   std::ostringstream preciceTraceStream; \
   preciceTraceStream << #var1 << "=" << var1; \
   precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define tpreciceTrace2(methodname, var1, var2) \
  std::string preciceMethodName(methodname); \
  std::ostringstream preciceTraceStream; \
  preciceTraceStream << #var1 << "=" << var1 << ", " << #var2 << "=" << var2; \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define tpreciceTrace3(methodname, var1, var2, var3) \
  std::string preciceMethodName(methodname); \
  std::ostringstream preciceTraceStream; \
  preciceTraceStream         << #var1 << "=" << var1  \
                     << ", " << #var2 << "=" << var2  \
                     << ", " << #var3 << "=" << var3; \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define tpreciceTrace4(methodname, var1, var2, var3, var4) \
  std::string preciceMethodName(methodname); \
  std::ostringstream preciceTraceStream; \
  preciceTraceStream         << #var1 << "=" << var1  \
                     << ", " << #var2 << "=" << var2  \
                     << ", " << #var3 << "=" << var3  \
                     << ", " << #var4 << "=" << var4; \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define tpreciceTrace5(methodname, var1, var2, var3, var4, var5) \
  std::string preciceMethodName(methodname); \
  std::ostringstream preciceTraceStream; \
  preciceTraceStream         << #var1 << "=" << var1  \
                     << ", " << #var2 << "=" << var2  \
                     << ", " << #var3 << "=" << var3  \
                     << ", " << #var4 << "=" << var4  \
                     << ", " << #var5 << "=" << var5; \
  precice::utils::Tracer preciceTracer(PRECICE_LOGGING_DEVICE, methodname, preciceTraceStream.str());

#define tpreciceTrace6(methodname, var1, var2, var3, var4, var5, var6) \
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

//#define tpreciceDebug(methodname, message)
#define tpreciceTrace(methodname)
#define tpreciceTrace1(methodname, var1)
#define tpreciceTrace2(methodname, var1, var2)
#define tpreciceTrace3(methodname, var1, var2, var3)
#define tpreciceTrace4(methodname, var1, var2, var3, var4)
#define tpreciceTrace5(methodname, var1, var2, var3, var4, var5)
#define tpreciceTrace6(methodname, var1, var2, var3, var4, var5, var6)

#endif // ! Debug

/**
 * @brief Needed for macros exiting program execution, to empty logging cache.
 */
//#ifdef _UTILS_LOG_H_
//#define PRECICE_CLOSE_LOGGER logging::Logger::criticalAbort();
//#else
//#define PRECICE_CLOSE_LOGGER
//#endif

/**
 * @brief Wrapper for logging::Logger error output with program exit.
 */
#define tpreciceError(methodname, message) \
   { \
      std::ostringstream conv; \
      conv << PRECICE_PROCESS_RANK_STREAM << " [PRECICE] ERROR: " << message; \
      PRECICE_LOGGING_DEVICE.error (methodname, conv.str()); \
      std::abort(); \
   }

#define tpreciceCheck(check, methodname, errormessage) \
   if ( !(check) ) { \
      tpreciceError(methodname, errormessage); \
   }


