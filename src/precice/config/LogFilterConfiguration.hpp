#ifndef _TARCH_LOGGING_CONFIGURATION_LOGFILTERCONFIGURATION_H_
#define _TARCH_LOGGING_CONFIGURATION_LOGFILTERCONFIGURATION_H_

#include "tarch/logging/Log.h"
#include "tarch/logging/CommandLineLogger.h"
#include "utils/xml/XMLTag.hpp"
#include <string>

namespace precice {
namespace config {

/**
 * @brief Configuration class for log-filters.
 */
class LogFilterConfiguration : public utils::XMLTag::Listener
{
public:

  LogFilterConfiguration ( utils::XMLTag& parent );

  //bool isValid() const;

  tarch::logging::CommandLineLogger::FilterList getFilterList() const;

  virtual void xmlTagCallback ( utils::XMLTag& tag );

  virtual void xmlEndTagCallback ( utils::XMLTag& tag );

private:

  static tarch::logging::Log _log;

  const std::string TAG;
  const std::string ATTR_TARGET;
  const std::string ATTR_COMPONENT;
  const std::string ATTR_SWITCH;

  //bool _isValid;

  tarch::logging::CommandLineLogger::FilterList _filterList;
};

}} // namespace precice, config

#endif
