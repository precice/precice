#include <boost/test/detail/log_level.hpp>
#include <boost/test/tools/fpc_tolerance.hpp>
#include <boost/test/tree/test_unit.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <iostream>

#include "logging/LogConfiguration.hpp"
#include "utils/ArgumentFormatter.hpp"

namespace precice::testing {

/// Boost test Initialization function
/**
Boost Test Log Levels and corresponding command line arguments to --log_level
as of Boost 1.68:

\code
type = enum boost::unit_test::log_level : int {
  boost::unit_test::invalid_log_level = -1,
  boost::unit_test::log_successful_tests = 0, all
  boost::unit_test::log_test_units = 1, unit_scope, test_suite
  boost::unit_test::log_messages = 2, message
  boost::unit_test::log_warnings = 3, warning
  boost::unit_test::log_all_errors = 4, error (default log level)
  boost::unit_test::log_cpp_exception_errors = 5, cpp_exception
  boost::unit_test::log_system_errors = 6, system_error
  boost::unit_test::log_fatal_errors, fatal_error
  boost::unit_test::log_nothing, nothing
}
\endcode
**/
boost::unit_test::log_level getBoostTestLogLevel()
{
  namespace bu = boost::unit_test;
#if BOOST_VERSION == 106900 || __APPLE__ && __MACH__
  std::cerr << "Boost 1.69 and macOS get log_level is broken, preCICE log level set to debug.\n";
  return bu::log_successful_tests;
#else
  return bu::runtime_config::get<bu::log_level>(bu::runtime_config::btrt_log_level);
#endif
}

std::string filterFromLogLevel(boost::unit_test::log_level logLevel)
{
  namespace bu = boost::unit_test;

  if (logLevel == bu::log_successful_tests || logLevel == bu::log_test_units) {
    return "";
  }
  if (logLevel == bu::log_messages) {
    return "%Severity% >= info";
  }
  if (logLevel == bu::log_warnings) {
    return "%Severity% >= warning";
  }
  // log warnings in any case
  return "%Severity% >= warning";
}

void setupTestLogging()
{
  // See if there is a manual override using a log.conf file.
  auto logConfigs = logging::readLogConfFile("log.conf");

  if (logConfigs.empty()) {
    // Nothing has been read from log.conf
    // We configure the log level based on the test framework log level

    logging::BackendConfiguration config;
    config.filter = filterFromLogLevel(getBoostTestLogLevel());

    const std::string prefix{"%TimeStamp(format=\"%H:%M:%S.%f\")%|%Participant%|%Rank%|%Module%|l%Line%|%Function%|"};

    // Console output
    config.format = prefix + "%ColorizedSeverity%%Message%";
    config.type   = "stream";
    config.output = "stdout";
    logConfigs.push_back(config);

    // File Outputs
    config.format = prefix + "%Severity%%Message%";
    config.type   = "file";
    config.output = "test.log";
    logConfigs.push_back(config);

    // The full debug log
    config.output = "test.debug.log";
    config.filter = "";
    logConfigs.push_back(config);
  }

  logging::setupLogging(logConfigs);
  // Lock the logging configuration to prevent systemtests from changing them
  logging::lockConf();
}

} // namespace precice::testing

// Fixtures need to be defined in the global scope

struct PreciceTestLoggingFixture {
  static void setup()
  {
    std::cerr << "Setup up logging\n";
    precice::testing::setupTestLogging();
  }
};

BOOST_TEST_GLOBAL_FIXTURE(PreciceTestLoggingFixture);
