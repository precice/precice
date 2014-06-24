// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#if !defined(_TARCH_LOGGING_CCA_LOGGER_H_) && defined(CCA)
#define _TARCH_LOGGING_CCA_LOGGER_H_


#include "tarch/multicore/BooleanSemaphore.h"


#include <string>


namespace tarch {
  namespace logging {
    class CCALogger;
  }
}



/**
 * CCA Logger
 *
 * Logger for ISCADT. It is basically a wrapper around the command line logger.
 *
 * @author Tobias Weinzierl
 */
class tarch::logging::CCALogger {
  public:
    /**
     * Port name for CCA.
     *
     * The time stepper can write data to a uses port. For this, the uses port
     * has to have this name, and you have to implement a UsesPortService that
     * can answer to this name. Usually your runner is this uses port service.
     */
    static const std::string CCAUsesPortLog;

  private:
    /**
     * Constructor
     *
     * The command line logger is a singleton. For the CCALogger, there's no
     * need to make it a singleton - the port stuff is a singleton anyway.
     * However, I made also the CCALogger a singleton to have the same
     * signature for the command line and the CCA stuff.
     */
    CCALogger();

    tarch::multicore::BooleanSemaphore _semaphore;

    /**
     * I'm having issues with recursive calls of the logging. Someone triggers
     * a log, the log tries to get a uses port, the uses port in turn tries to
     * log, and so forth. We can avoid this as the cca logger is a singleton.
     * This flag is set whenever a function is called. If we detect a recursive
     * call, we do not try to forward data to the uses port, but we write the
     * output directly to the terminal.
     */
    bool _isCalled;
  public:
    ~CCALogger();

    static CCALogger& getInstance();

    void debug(      const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message);
    void info(       const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message);
    void warning(    const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message);
    void error(      const long int& timestampMS, const std::string& timestampHumanReadable, const std::string& machineName, const std::string& trace, const std::string& message);

    void indent( bool indent );
};

#endif
