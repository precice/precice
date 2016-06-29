#pragma once

#include "logging/LogConfiguration.hpp"
#include "utils/xml/XMLTag.hpp"

namespace precice {
namespace config {

/// Configures the log config file to use
class LogConfiguration : public utils::XMLTag::Listener
{
public:
  LogConfiguration(utils::XMLTag& parent);

  virtual void xmlTagCallback(utils::XMLTag& tag);

  virtual void xmlEndTagCallback(utils::XMLTag& tag);

private:
  static precice::logging::Logger _log;

  precice::logging::LoggingConfiguration _logconfig;
};

}} // namespace precice, config
