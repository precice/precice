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
  std::string type = "stream";
  std::string output = "stdout";
  boost::log::filter filter = boost::log::parse_filter("%Severity% > debug");
  boost::log::basic_formatter<char> format = boost::log::parse_formatter("(%Rank%) %TimeStamp(format=\"%H:%M:%S\")% [%Module%]:%Line% in %Function%: %Message%");

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
