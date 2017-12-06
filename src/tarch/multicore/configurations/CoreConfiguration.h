#ifndef _TARCH_MULTICORE_TBB_CONFIGURATION_CORE_CONFIGURATION_H_
#define _TARCH_MULTICORE_TBB_CONFIGURATION_CORE_CONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "logging/Logger.hpp"
#include <string>


namespace tarch {
  namespace multicore {
    namespace configurations {
      class CoreConfiguration;
    }
  }
}


/**
 * Core Configuration
 *
 * @author Tobias Weinzierl
 */
class tarch::multicore::configurations::CoreConfiguration: public tarch::configuration::Configuration {
  private:
    /**
     * Log device.
     */
    static precice::logging::Logger _log;

    int _numberOfThreads;

    bool _hasParsed;

  public:
    CoreConfiguration();
    virtual ~CoreConfiguration();

    virtual std::string getTag() const;

	virtual void parseSubtag( precice::xml::Parser::CTag *pTag);

    virtual bool isValid() const;

    virtual void toXML(std::ostream& out) const;

    /**
     * @return Number of threads or 0 if you should use the standard number.
     */
    int getNumberOfThreads() const;
};


#endif
