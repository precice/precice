#pragma once

#include "logging/LogConfiguration.hpp"
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice::logging {

/// Configures the log config file to use
class LogConfiguration : public xml::XMLTag::Listener {
public:
  LogConfiguration(xml::XMLTag &parent);

  void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag) override;

  void xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &tag) override;

private:
  precice::logging::Logger _log{"logging::config::LogConfiguration"};

  precice::logging::LoggingConfiguration _logconfig;
};

} // namespace precice::logging
