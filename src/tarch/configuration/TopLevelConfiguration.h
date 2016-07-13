#ifndef _TARCH_CONFIGURATION_TOP_LEVEL_CONFIGURATION_H_
#define _TARCH_CONFIGURATION_TOP_LEVEL_CONFIGURATION_H_

#ifdef Parallel
#include <mpi.h>
#endif

#include "tarch/configuration/Configuration.h"


#define registerTopLevelConfiguration(name) \
  static tarch::configuration::TopLevelConfigurationFactory<name> thisConfigurationInstanceFactory;


namespace tarch {
  namespace configuration {
    class TopLevelConfiguration;
  }
}

/**
 * Abstract supertype of all interfaces.
 *
 * !!! Rationale
 *
 * - Do not use a logger within a top level configuration. Otherwise, the code
 *   might crash.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.21 $
 */
class tarch::configuration::TopLevelConfiguration: public tarch::configuration::Configuration {
  public:
    /**
     * Destructor
     */
    virtual ~TopLevelConfiguration();

    /**
     * Clone configuration. Receiver is responsible to destroy the instance himself.
     */
    virtual TopLevelConfiguration* clone() const = 0;

    /**
     * Take this configuration and run it. If the configuration returns 0,
     * everything is fine. A value greater than 0 is an error code and makes
     * the application terminate.
     */
    virtual int interpreteConfiguration() = 0;
};


#endif

