// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LOGGING_CONFIGURATION_LOGFILTERCONFIGURATION_H_
#define _TARCH_LOGGING_CONFIGURATION_LOGFILTERCONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/logging/Log.h"
#include "tarch/logging/CommandLineLogger.h"
#include <string>

namespace tarch {
  namespace logging {
    namespace configurations {
      class LogFilterConfiguration;
    }
  }
}


/**
 * @author Tobias Weinzierl
 */
class tarch::logging::configurations::LogFilterConfiguration: public tarch::configuration::Configuration {
  private:
    /**
     * Log device.
     */
    static tarch::logging::Log _log;

    /**
     * Represents one entry of the filter list. Syntax:
     *
     * <log-filter target="debug" component="component-name" rank="int|*" switch="on|off" />
     */
    tarch::logging::CommandLineLogger::FilterList _filterList;

    /**
     * Holds true if the filter configuration is valid.
     */
    bool _isValid;

  public:
    LogFilterConfiguration();
    virtual ~LogFilterConfiguration();

    tarch::logging::CommandLineLogger::FilterList getFilterList() const;

    virtual std::string getTag() const;

    virtual void parseSubtag( tarch::irr::io::IrrXMLReader* _xmlReader );

    virtual bool isValid() const;

    virtual void toXML(std::ostream& out) const;
};

#endif
