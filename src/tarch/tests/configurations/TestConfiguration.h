// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_TESTS_CONFIGURATION_TESTCONFIGURATION_H_
#define _TARCH_TESTS_CONFIGURATION_TESTCONFIGURATION_H_

#include "tarch/configuration/TopLevelConfiguration.h"
#include "tarch/logging/Log.h"
#include "tarch/logging/configurations/LogOutputFormatConfiguration.h"
#include "tarch/logging/configurations/LogFilterConfiguration.h"
#include <string>

namespace tarch {
  namespace tests {
    namespace configurations {
      class TestConfiguration;
    }
  }
}


/**
 * There's no opportunity to configure the automatic tests, but for consistency
 * reasons there's a configuration class for this use case. For every instance
 * of this class the test suite has to be run.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.2 $
 */
class tarch::tests::configurations::TestConfiguration: public tarch::configuration::TopLevelConfiguration {
  private:
    /**
     * Log interface the class writes to.
     */
    static tarch::logging::Log _log;

    bool _isValid;

    tarch::logging::configurations::LogFilterConfiguration _logConfiguration;
    tarch::logging::configurations::LogOutputFormatConfiguration _logFormatConfiguration;
  public:
    /**
     * Constructor.
     */
    TestConfiguration();

    /**
     * Destructor.
     */
    virtual ~TestConfiguration();

    virtual std::string getTag() const;

    virtual void parseSubtag( tarch::irr::io::IrrXMLReader* _xmlReader );

    virtual bool isValid() const;

    virtual void toXML(std::ostream& out) const;

    virtual tarch::configuration::TopLevelConfiguration* clone() const;

    virtual int interpreteConfiguration();
};

#endif
