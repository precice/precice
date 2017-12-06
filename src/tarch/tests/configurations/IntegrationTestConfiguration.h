#ifndef _TARCH_TESTS_CONFIGURATION_INTEGRATION_TESTCONFIGURATION_H_
#define _TARCH_TESTS_CONFIGURATION_INTEGRATION_TESTCONFIGURATION_H_

#include "tarch/configuration/TopLevelConfiguration.h"
#include "logging/Logger.hpp"
#include "tarch/logging/configurations/LogFilterConfiguration.h"
#include "tarch/logging/configurations/LogOutputFormatConfiguration.h"
#include <string>

namespace tarch {
  namespace tests {
    namespace configurations {
      class IntegrationTestConfiguration;
    }
  }
}


/**
 * @author Tobias Weinzierl
 * @version $Revision: 1.2 $
 */
class tarch::tests::configurations::IntegrationTestConfiguration: public tarch::configuration::TopLevelConfiguration {
  private:
    /**
     * Log interface the class writes to.
     */
    static precice::logging::Logger _log;

    bool _isValid;

    std::string _outputDirectoryForTempFiles;

    tarch::logging::configurations::LogFilterConfiguration       _logConfiguration;
    tarch::logging::configurations::LogOutputFormatConfiguration _logFormatConfiguration;
  public:
    /**
     * Constructor.
     */
    IntegrationTestConfiguration();

    /**
     * Destructor.
     */
    virtual ~IntegrationTestConfiguration();

    virtual std::string getTag() const;

	virtual void parseSubtag( precice::xml::Parser::CTag *pTag);
	
    virtual bool isValid() const;

    virtual void toXML(std::ostream& out) const;

    virtual TopLevelConfiguration* clone() const;

    virtual int interpreteConfiguration();
};

#endif
