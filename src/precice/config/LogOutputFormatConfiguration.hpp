// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CONFIG_LOGOUTPUTFORMATCONFIGURATION_HPP
#define PRECICE_CONFIG_LOGOUTPUTFORMATCONFIGURATION_HPP

#include "utils/xml/XMLTag.hpp"
#include "tarch/logging/Log.h"
#include <string>

namespace precice {
namespace config {

/**
 * @brief Configuration for log output format.
 */
class LogOutputFormatConfiguration : public utils::XMLTag::Listener
{
public:

  LogOutputFormatConfiguration ( utils::XMLTag& parent );

  //const std::string& getTag() const;

  //virtual void parseSubtag ( tarch::irr::io::IrrXMLReader* _xmlReader );

  //bool isValid() const;

  const std::string& getLogColumnSeparator() const;

  bool getLogTimeStamp() const;

  bool getLogTimeStampHumanReadable() const;

  bool getLogMachineName() const;

  bool getLogMessageType() const;

  bool getLogTrace() const;

  virtual void xmlTagCallback ( utils::XMLTag& tag );

  virtual void xmlEndTagCallback ( utils::XMLTag& tag );

private:

  static tarch::logging::Log _log;

  const std::string TAG;
  const std::string ATTR_COLUMN_SEPARATOR;
  const std::string ATTR_LOG_TIMESTAMP;
  const std::string ATTR_LOG_TIMESTAMP_HUMAN_READABLE;
  const std::string ATTR_LOG_MACHINE_NAME;
  const std::string ATTR_LOG_MESSAGE_TYPE;
  const std::string ATTR_LOG_TRACE;

  std::string _logColumnSeparator;

  bool _logTimeStamp;

  bool _logTimeStampHumanReadable;

  bool _logMachineName;

  bool _logMessageType;

  bool _logTrace;

  //bool _isValid;
};

}} // namespace precice, config

#endif // PRECICE_CONFIG_LOGOUTPUTFORMATCONFIGURATION_HPP
