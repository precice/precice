#ifndef _TARCH_LOGGING_CONFIGURATION_LOG_OUTPUT_FORMAT_CONFIGURATION_H_
#define _TARCH_LOGGING_CONFIGURATION_LOG_OUTPUT_FORMAT_CONFIGURATION_H_

#include "tarch/configuration/Configuration.h"
#include "tarch/logging/Log.h"
#include <string>

namespace tarch {
  namespace logging {
    namespace configurations {
      class LogOutputFormatConfiguration;
    }
  }
}


/**
 * @author Tobias Weinzierl
 */
class tarch::logging::configurations::LogOutputFormatConfiguration: public tarch::configuration::Configuration {
  private:
    /**
     * Log device.
     */
    static tarch::logging::Log _log;

    std::string  _logColumnSeparator;
    bool         _logTimeStamp;
    bool         _logTimeStampHumanReadable;
    bool         _logMachineName;
    bool         _logMessageType;
    bool         _logTrace;
    std::string  _logOutputFile;
    bool         _isValid;
    bool         _hasParsed;
  public:
    LogOutputFormatConfiguration();
    virtual ~LogOutputFormatConfiguration();

    virtual std::string getTag() const;

	virtual void parseSubtag( precice::xml::ConfigParser::CTag *pTag);
	
    virtual bool isValid() const;

    virtual void toXML(std::ostream& out) const;

    bool hasParsed() const;

    std::string getLogColumnSeparator() const;
    bool getLogTimeStamp() const;
    bool getLogTimeStampHumanReadable() const;
    bool getLogMachineName() const;
    bool getLogMessageType() const;
    bool getLogTrace() const;
    std::string getLogFileName() const;
};

#endif
