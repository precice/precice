#pragma once

#include <string>
#include <vector>
#include <boost/log/utility/setup/filter_parser.hpp>
#include <boost/log/utility/setup/formatter_parser.hpp>

namespace precice {
namespace logging {


/// Holds the configuration for one logging backend (sink) and takes care of default values.
struct BackendConfiguration
{
  static const std::string default_type;
  static const std::string default_output;
  static const std::string default_filter;
  static const std::string default_formatter;
    
  std::string type = default_type;
  std::string output = default_output;
  std::string filter = default_filter;
  std::string format = default_formatter;
  bool enabled = true;

  /// Sets on option, overwrites default values.
  void setOption(std::string key, std::string value);
};

/// Holds the configuration of the logging system
using LoggingConfiguration = std::vector<BackendConfiguration>;

void setupLogging(std::string logConfigFile = "log.conf");
void setupLogging(LoggingConfiguration configs);

/// Sets the current MPI rank as a logging attribute
void setMPIRank(const int rank);

}} // namespace precice, logging
