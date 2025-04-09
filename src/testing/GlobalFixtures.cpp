#include <boost/test/tree/test_unit.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <iostream>

#include "logging/LogConfiguration.hpp"
#include "utils/ArgumentFormatter.hpp"

namespace precice::testing {

void setupTestLogging()
{
  // See if there is a manual override using a log.conf file.
  auto logConfigs = logging::readLogConfFile("log.conf");

  if (logConfigs.empty()) {
    // Nothing has been read from log.conf
    // We configure the log level based on the test framework log level

    logging::BackendConfiguration config;
    config.filter = "%Severity% >= debug";

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
    config.filter = "%Severity% >= debug";
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
