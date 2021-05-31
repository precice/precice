#pragma once

#include <string>
#include <vector>

namespace precice {
namespace logging {

/// Holds the configuration for one logging backend (sink) and takes care of default values.
struct BackendConfiguration {
  static const std::string default_type;
  static const std::string default_output;
  static const std::string default_filter;
  static const std::string default_formatter;

  std::string type    = default_type;
  std::string output  = default_output;
  std::string filter  = default_filter;
  std::string format  = default_formatter;
  bool        enabled = true;

  /// Sets on option, overwrites default values.
  void setOption(std::string key, std::string value);
};

/// Holds the configuration of the logging system
using LoggingConfiguration = std::vector<BackendConfiguration>;

/// Reads a log configuration file, returns vector of BackEndConfiguration
LoggingConfiguration readLogConfFile(std::string const &filename);

/// Configures the logging from a log file
void setupLogging(std::string const &logConfigFile = "log.conf");

/// Configures the logging from a LoggingConfiguration
void setupLogging(LoggingConfiguration configs, bool enabled = true);

/// Sets the current MPI rank as a logging attribute
void setMPIRank(int const rank);

/// Sets the name of the current participant as a logging attribute
void setParticipant(std::string const &name);

/// Locks the configuration, ignoring any future calls to setupLogging()
void lockConf();

/** The global lock for the log configuration.
 * This variable is always initialized to false and can only be modified by calling lockConf()
 */
extern bool _precice_logging_config_lock;

} // namespace logging
} // namespace precice
