#pragma once

#include "logging/LogConfiguration.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace config {

/// Configures the log config file to use
class LogConfiguration : public xml::XMLTag::Listener
{
public:
  LogConfiguration(xml::XMLTag& parent);

  virtual void xmlTagCallback(xml::XMLTag& tag);

  virtual void xmlEndTagCallback(xml::XMLTag& tag);

private:
  static precice::logging::Logger _log;

  precice::logging::LoggingConfiguration _logconfig;
};

}} // namespace precice, config
