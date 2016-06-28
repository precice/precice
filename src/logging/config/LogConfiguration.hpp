#pragma once

#include "logging/Logger.hpp"
#include "utils/xml/XMLTag.hpp"
#include <string>

namespace precice {
namespace config {

/// Configures the log config file to use
class LogConfiguration : public utils::XMLTag::Listener
{
public:
  LogConfiguration(utils::XMLTag& parent);

  virtual void xmlTagCallback ( utils::XMLTag& tag );

  virtual void xmlEndTagCallback ( utils::XMLTag& tag );

private:
  static precice::logging::Logger _log;

  const std::string TAG;
  const std::string ATTR_FILE;
};

}} // namespace precice, config
