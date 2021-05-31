#pragma once

#include <string>
#include "logging/LogConfiguration.hpp"
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace config {

/// Configures the log config file to use
class LogConfiguration : public xml::XMLTag::Listener {
public:
  LogConfiguration(xml::XMLTag &parent);

  virtual void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag);

  virtual void xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag);

private:
  precice::logging::Logger _log{"logging::config::LogConfiguration"};

  precice::logging::LoggingConfiguration _logconfig;
};

} // namespace config
} // namespace precice
