// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PARALLEL_CONFIGURATION_PARALLELCONFIGURATION_H_
#define _TARCH_PARALLEL_CONFIGURATION_PARALLELCONFIGURATION_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include "tarch/logging/Log.h"
#include "tarch/configuration/Configuration.h"

#ifdef Parallel
#endif

namespace tarch {
  namespace parallel {
    namespace configuration {
      class ParallelConfiguration;
    }
  }
}


/**
 * Parallel configuration
 *
 * The parallel configuration holds all the attributes for configuring the
 * parallel component's classes, i.e. it should be hold by any configuration
 * that identifies a parallel runner. It configures the deadlocks, timeout
 * warning, communicators, and so forth. The
 * validation returns an error if no value is set and -DParallel is set.
 * Furthermore, the operation yields a warning if -DParallel is not set but a
 * parallel tag is configured.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.19 $
 */
class tarch::parallel::configuration::ParallelConfiguration: public tarch::configuration::Configuration {
  private:
    /**
     * Log device.
     */
    static tarch::logging::Log _log;

    /**
     * Specifies the period until a waiting process writes a warning (in
     * seconds).
     */
    clock_t _timeoutWarning;

    /**
     * Specifies the period until a waiting process shuts down (in seconds).
     */
    clock_t _deadlockTimeOut;

    /**
     * 'default' means take the standard communicator. This value of the xml
     * file is represented by -1.
     */
    #ifdef Parallel
    MPI_Comm _communicator;
    #endif

    bool _validCommunicator;

    static const std::string WaitWarningTime;
    static const std::string DeadlockTimeout;
    static const std::string Communicator;

    bool writeWarning() const;
  public:
  	ParallelConfiguration();

	  virtual ~ParallelConfiguration();


    virtual std::string getTag() const;

    /**
     * Parse a TAG-tag. The parser reads the tag, validates the context, ensures
     * the corresponding compiler switch is set or not (depending on tag) and
     * returns. The argument may not be 0. If the either the validation or the
     * compiler switch check fails, all successing isValid() calls fail.
     */
    void parseSubtag( irr::io::IrrXMLReader* _xmlReader );

    /**
     * @returns Whether all compiler switches checked using the parseSubtag()
     *          routine were valid compared to the values specified in the
     *          configuration file.
     */
    bool isValid() const;

    virtual void toXML(std::ostream& out) const;

    void interpreteConfiguration() const;
};

#endif
