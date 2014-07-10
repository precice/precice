// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LOGGING_LOG_H_
#define _TARCH_LOGGING_LOG_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include <iostream>
#include <sstream>


namespace tarch {
  namespace logging {
    class Log;
    class CommandLineLogger;
    class CCALogger;
  }
}

#ifdef Debug
/**
 * @see logInfo() macro
 */
#define logDebug(methodName, messageStream) \
   { \
      std::ostringstream conv; \
      conv << messageStream; \
      conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
      _log.debug (methodName, conv.str()); \
   }

#define logDebugMasterOnly(methodName, messageStream) \
   { \
      std::ostringstream conv; \
      conv << messageStream; \
      conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
      _log.debugMasterOnly (methodName, conv.str()); \
   }

#define logTraceIn(methodName) \
  { \
    std::ostringstream conv; \
    conv << "in (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceInWith1Argument(methodName,argument0) \
  { \
    std::ostringstream conv; \
    conv << "in:" << #argument0 << ":" << argument0; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceInWith2Arguments(methodName,argument0,argument1) \
  { \
    std::ostringstream conv; \
    conv << "in:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceInWith3Arguments(methodName,argument0,argument1,argument2) \
  { \
    std::ostringstream conv; \
    conv << "in:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceInWith4Arguments(methodName,argument0,argument1,argument2,argument3) \
  { \
    std::ostringstream conv; \
    conv << "in:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceInWith5Arguments(methodName,argument0,argument1,argument2,argument3,argument4) \
  { \
    std::ostringstream conv; \
    conv << "in:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << "," << #argument4 << ":" << argument4; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceInWith6Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5) \
  { \
    std::ostringstream conv; \
    conv << "in:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << "," << #argument4 << ":" << argument4; \
    conv << "," << #argument5 << ":" << argument5; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceInWith7Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6) \
  { \
    std::ostringstream conv; \
    conv << "in:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << "," << #argument4 << ":" << argument4; \
    conv << "," << #argument5 << ":" << argument5; \
    conv << "," << #argument6 << ":" << argument6; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceInWith8Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6,argument7) \
  { \
    std::ostringstream conv; \
    conv << "in:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << "," << #argument4 << ":" << argument4; \
    conv << "," << #argument5 << ":" << argument5; \
    conv << "," << #argument6 << ":" << argument6; \
    conv << "," << #argument7 << ":" << argument7; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceInWith9Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6,argument7,argument8) \
  { \
    std::ostringstream conv; \
    conv << "in:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << "," << #argument4 << ":" << argument4; \
    conv << "," << #argument5 << ":" << argument5; \
    conv << "," << #argument6 << ":" << argument6; \
    conv << "," << #argument7 << ":" << argument7; \
    conv << "," << #argument8 << ":" << argument8; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.debug (methodName, conv.str()); \
    _log.indent(true,_log.getTraceInformation(methodName),conv.str()); \
  }

#define logTraceOut(methodName) \
  { \
    std::ostringstream conv; \
    conv << "out (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }

#define logTraceOutWith1Argument(methodName,argument0) \
  { \
    std::ostringstream conv; \
    conv << "out:" << #argument0 << ":" << argument0; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }

#define logTraceOutWith2Arguments(methodName,argument0,argument1) \
  { \
    std::ostringstream conv; \
    conv << "out:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }

#define logTraceOutWith3Arguments(methodName,argument0,argument1,argument2) \
  { \
    std::ostringstream conv; \
    conv << "out:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }

#define logTraceOutWith4Arguments(methodName,argument0,argument1,argument2,argument3) \
  { \
    std::ostringstream conv; \
    conv << "out:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }

#define logTraceOutWith5Arguments(methodName,argument0,argument1,argument2,argument3,argument4) \
  { \
    std::ostringstream conv; \
    conv << "out:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << "," << #argument4 << ":" << argument4; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }

#define logTraceOutWith6Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5) \
  { \
    std::ostringstream conv; \
    conv << "out:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << "," << #argument4 << ":" << argument4; \
    conv << "," << #argument5 << ":" << argument5; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }

#define logTraceOutWith7Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6) \
  { \
    std::ostringstream conv; \
    conv << "out:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << "," << #argument4 << ":" << argument4; \
    conv << "," << #argument5 << ":" << argument5; \
    conv << "," << #argument6 << ":" << argument6; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }


#define logTraceOutWith8Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6,argument7) \
  { \
    std::ostringstream conv; \
    conv << "out:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1 << ":" << argument1; \
    conv << "," << #argument2 << ":" << argument2; \
    conv << "," << #argument3 << ":" << argument3; \
    conv << "," << #argument4 << ":" << argument4; \
    conv << "," << #argument5 << ":" << argument5; \
    conv << "," << #argument6 << ":" << argument6; \
    conv << "," << #argument7 << ":" << argument7; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }

#define logTraceOutWith12Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6,argument7,argument8,argument9,argument10,argument11) \
  { \
    std::ostringstream conv; \
    conv << "out:" << #argument0 << ":" << argument0; \
    conv << "," << #argument1  << ":" << argument1; \
    conv << "," << #argument2  << ":" << argument2; \
    conv << "," << #argument3  << ":" << argument3; \
    conv << "," << #argument4  << ":" << argument4; \
    conv << "," << #argument5  << ":" << argument5; \
    conv << "," << #argument6  << ":" << argument6; \
    conv << "," << #argument7  << ":" << argument7; \
    conv << "," << #argument8  << ":" << argument8; \
    conv << "," << #argument9  << ":" << argument9; \
    conv << "," << #argument10 << ":" << argument10; \
    conv << "," << #argument11 << ":" << argument11; \
    conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
    _log.indent(false,_log.getTraceInformation(methodName),conv.str()); \
    _log.debug (methodName, conv.str()); \
  }

#else
#define logDebug(methodName, messageStream)
#define logDebugMasterOnly(methodName, messageStream)
#define logTraceIn(methodName)
#define logTraceInWith1Argument(methodName,argument0)
#define logTraceInWith2Arguments(methodName,argument0,argument1)
#define logTraceInWith3Arguments(methodName,argument0,argument1,argument2)
#define logTraceInWith4Arguments(methodName,argument0,argument1,argument2,argument3)
#define logTraceInWith5Arguments(methodName,argument0,argument1,argument2,argument3,argument4)
#define logTraceInWith6Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5)
#define logTraceInWith7Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6)
#define logTraceInWith8Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6,argument7)
#define logTraceInWith9Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6,argument7,argument8)
#define logTraceOut(methodName)
#define logTraceOutWith1Argument(methodName,argument0)
#define logTraceOutWith2Arguments(methodName,argument0,argument1)
#define logTraceOutWith3Arguments(methodName,argument0,argument1,argument2)
#define logTraceOutWith4Arguments(methodName,argument0,argument1,argument2,argument3)
#define logTraceOutWith5Arguments(methodName,argument0,argument1,argument2,argument3,argument4)
#define logTraceOutWith6Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5)
#define logTraceOutWith7Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6)
#define logTraceOutWith8Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6,argument7)
#define logTraceOutWith12Arguments(methodName,argument0,argument1,argument2,argument3,argument4,argument5,argument6,argument7,argument8,argument9,argument10,argument11)
#endif


/**
 * @brief Wrapper macro around tarch::tarch::logging::Log to improve logging.
 *
 * A Log object with name _log has to be defined at the place of calling this
 * macro.
 *
 * This macro allows to combine strings and variables arbitrarily
 * in an efficient way (only one ostringstream object has to be created per
 * usage of logInfo).
 *
 * Usage:
 * logInfo( "myOperation()", "anyText" << myVar << ",anotherText" << myVar2 );
 *
 * !!! Hint
 *
 * Never use the + operator to concatenate data as this error is error-prone.
 * If you use always the << operator, you are on the save side, as the +
 * operator works only for strings properly. If you use it with a string and
 * another data type, it might be that the string is assigned an invalid length.
 */
#define logInfo(methodName, messageStream) \
   { \
      std::ostringstream conv; \
      conv << messageStream; \
      _log.info (methodName, conv.str()); \
   }

#define logInfoMasterOnly(methodName, messageStream) \
   { \
      std::ostringstream conv; \
      conv << messageStream; \
      _log.infoMasterOnly (methodName, conv.str()); \
   }

#define logExceptionAndQuit(exception) \
  { \
    std::cerr << std::string("caught exception (file:)") << __FILE__ << std::string(",line:") << __LINE__ << std::string("): ") << std::string(exception.what()); \
    exit(-1); \
  }


/**
 * @brief Wrapper macro around tarch::tarch::logging::Log to improve logging.
 *
 * A Log object with name _log has to be defined at the place of calling this
 * macro.
 *
 * This macro allows to combine strings and variables arbitrarily
 * in an efficient way (only one ostringstream object has to be created per
 * usage of logWarning).
 *
 * Usage:
 * logWarning( "myOperation()", "anyText" << myVar << ",anotherText" << myVar2 );
 */
#define logWarning(methodName, messageStream) \
   { \
      std::ostringstream conv; \
      conv << messageStream; \
      conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
      _log.warning (methodName, conv.str()); \
   }

/**
 * @brief Wrapper macro around tarch::tarch::logging::Log to improve logging.
 *
 * A Log object with name _log has to be defined at the place of calling this
 * macro.
 *
 * This macro allows to combine strings and variables arbitrarily
 * in an efficient way (only one ostringstream object has to be created per
 * usage of logWarningMasterOnly).
 *
 * Usage:
 * logWarning( "myOperation()", "anyText" << myVar << ",anotherText" << myVar2 );
 */
#define logWarningMasterOnly(methodName, messageStream) \
   { \
      std::ostringstream conv; \
      conv << messageStream; \
      conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
      _log.warningMasterOnly (methodName, conv.str()); \
   }

/**
 * @brief Wrapper macro around tarch::tarch::logging::Log to improve logging.
 *
 * A Log object with name _log has to be defined at the place of calling this
 * macro.
 *
 * This macro allows to combine strings and variables arbitrarily
 * in an efficient way (only one ostringstream object has to be created per
 * usage of logError).
 *
 * Usage:
 * logInfo( "myOperation()", "anyText" << myVar << ",anotherText" << myVar2 );
 */
#define logError(methodName, messageStream) \
   { \
      std::ostringstream conv; \
      conv << messageStream; \
      conv << " (file:" << __FILE__ << ",line:" << __LINE__ << ")"; \
      _log.error (methodName, conv.str()); \
   }



/**
 * Log Device
 *
 * Log is the class all logging classes should use. To use the logging api they
 * have to create an instance by their own. It is suggested to hold this
 * instance static for the constructor of the Log class has to be given the
 * class name of the logging class. The Log class itself is stateless.
 * The log requests on this instance are processed here and forwarded to the
 * assigned logger (an internal attribute).
 *
 * Which concrete implementation has to be used for logging is switched using
 * a compiler attribute. Since the logging is used extremly often, this is
 * better than dynamic binding. Furthermoren, the log device is configured
 * using some compiler switches that define the type of information written.
 *
 * There are five different log levels, the user may write any output:
 *
 * - error: Here only errors have to be written to. Error messages may never
 *          be oppressed.
 * - warning:
 * - debug: Debug information that is switched off normally.
 * - info:  Statistical information, copyright and similar information. Should
 *          be used rather seldom.
 *
 * IMPORTANT: The underlying log device (like the CommandlineLogger) has to offer
 *            synchronized output methods. Thus, calls to logging methdos in
 *            multithreaded environments mark synchronization points of your programme
 *            and will change timing behaviour!
 *
 *
 * @version $Revision: 1.13 $
 * @author  Tobias Weinzierl
 */
class tarch::logging::Log {
  private:
    #ifdef CCA
    typedef CCALogger UsedLogService;
    #else
    typedef CommandLineLogger UsedLogService;
    #endif

    /**
     * Writes a timestamp to the standard output.
     */
    long int getTimeStampMS() const;

    std::string getTimeStampHumanReadable() const;

    /**
     * Name of the class that is using the interface.
     */
    std::string _className;
  public:
    /**
     * Writes information about the computer the output is written from.
     * The information string contains the operation system name, the computer
     * name and the cpu name.
     */
    std::string getMachineInformation() const;

    /**
     * Constructor
     *
     * @param className Name of the class that is using the logging component.
     *                  Please specify both class name and namespace using the
     *                  format namespace::classname
     */
    Log(const std::string& className);

    /**
     * Destructor
     */
    virtual ~Log();

    /**
     * Logs and exception and stops the application. I recommend to use the
     * corresponding macro instead.
     */
    static void exception( const std::bad_alloc& exception, const std::string& file, const int lineNumber );

    /**
     * Log Debug Information
     *
     * Remark: Please do not use this operation directly, but use the macros
     * above instead.
     *
     * The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    #if defined(Debug) && !defined(LogOff)
    void debug(const std::string& methodName, const std::string& message) const;
    #else
    void debug(const std::string& methodName, const std::string& message) const {
    }
    #endif

    /**
     * Log Debug Information
     *
     * Remark: Please do not use this operation directly, but use the macros
     * above instead.
     *
     * The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    #if defined(Debug) && !defined(LogOff)
    void debugMasterOnly(const std::string& methodName, const std::string& message) const;
    #else
    void debugMasterOnly(const std::string& methodName, const std::string& message) const {
    }
    #endif

    /**
     * Log Information
     *
     * The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    #if defined(LogOff)
    void info(const std::string& methodName, const std::string& message) const {}
    #else
    void info(const std::string& methodName, const std::string& message) const;
    #endif

    /**
     * Log Information
     *
     * The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    #if defined(LogOff)
    void infoMasterOnly(const std::string& methodName, const std::string& message) const {}
    #else
    void infoMasterOnly(const std::string& methodName, const std::string& message) const;
    #endif


    /**
     * Log a Warning
     *
     * The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    #if defined(LogOff)
    void warning(const std::string& methodName, const std::string& message) const {}
    #else
    void warning(const std::string& methodName, const std::string& message) const;
    #endif

    /**
     * Log a Warning
     *
     * The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    #if defined(LogOff)
    void warningMasterOnly(const std::string& methodName, const std::string& message) const {}
    #else
    void warningMasterOnly(const std::string& methodName, const std::string& message) const;
    #endif



    /**
     * Log an Error
     *
     * The method has to be given the method name. Often it
     * is useful to pass the whole method signature, that means e.g. if the
     * method info itself would log, the method invocation would look like this:
     * <pre>
     *   info("info(std::string,std::string)", "just an example text");
     * </pre>
     * So you are able to identify methods in your log although you use
     * overloading.
     *
     * @param methodName method name
     * @param message    log message
     */
    #if defined(LogOff)
    void error(const std::string& methodName, const std::string& message) const {}
    #else
    void error(const std::string& methodName, const std::string& message) const;
    #endif

    /**
     * Indent the Subsequent Messages
     *
     * Depending on indent the operation increments or decrements the indent.
     * The call is forwarded to the logging device.
     */
    #if defined(LogOff)
    void indent( bool indent, const std::string& trace, const std::string& message ) const {}
    #else
    void indent( bool indent, const std::string& trace, const std::string& message ) const;
    #endif

    std::string getTraceInformation( const std::string& methodName ) const;
};


#endif
