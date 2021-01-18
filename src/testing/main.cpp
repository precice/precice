#include <boost/test/tree/test_case_counter.hpp>
#include <boost/test/tree/traverse.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "com/SharedPointer.hpp"
#include "logging/LogConfiguration.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

namespace precice {
extern bool syncMode;
} // namespace precice

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
bool init_unit_test()
{
  using namespace boost::unit_test;
  using namespace precice;

  auto &master_suite        = framework::master_test_suite();
  master_suite.p_name.value = "preCICE Tests";

  auto logConfigs = logging::readLogConfFile("log.conf");

  if (logConfigs.empty()) { // nothing has been read from log.conf
#if BOOST_VERSION == 106900 || __APPLE__ && __MACH__
    std::cerr << "Boost 1.69 and macOS get log_level is broken, preCICE log level set to debug.\n";
    auto logLevel = log_successful_tests;
#else
    auto logLevel = runtime_config::get<log_level>(runtime_config::btrt_log_level);
#endif

    logging::BackendConfiguration config;
    if (logLevel == log_successful_tests or logLevel == log_test_units)
      config.filter = "%Severity% >= debug";
    if (logLevel == log_messages)
      config.filter = "%Severity% >= info";
    if (logLevel == log_warnings)
      config.filter = "%Severity% >= warning";
    if (logLevel >= log_all_errors)
      config.filter = "%Severity% >= warning"; // log warnings in any case

    const std::string prefix{"%TimeStamp(format=\"%H:%M:%S.%f\")%|%Participant%|%Rank%|%Module%|l%Line%|%Function%|"};

    // Console output
    config.format = prefix + "%ColorizedSeverity%%Message%";
    config.type   = "stream";
    config.output = "stdout";
    logConfigs.push_back(config);

    // File Outputs
    config.format = prefix + "%Severity%%Message%";
    config.type   = "file";

    // Same as console output
    config.output = "test.log";
    logConfigs.push_back(config);

    // The full debug log
    config.output = "test.debug.log";
    config.filter = "%Severity% >= debug";
    logConfigs.push_back(config);
  }

  // Initialize either with empty logConfig -> default, configs that are read from file
  // or from the Boost Test log level.
  logging::setupLogging(logConfigs);
  logging::lockConf();

  // Sets the default tolerance for floating point comparisions
  // Can be overwritten on a per-test or per-suite basis using decators
  // boost::unit_test::decorator::collector::instance() * boost::unit_test::tolerance(0.001);
  *tolerance(1e-9); // Stores the decorator in the collector singleton
#if BOOST_VERSION < 106900
  decorator::collector::instance().store_in(master_suite);
#else
  decorator::collector_t::instance().store_in(master_suite);
#endif

  return true;
}

int countEnabledTests()
{
  using namespace boost::unit_test;
  test_case_counter tcc;
  traverse_test_tree(framework::master_test_suite(), tcc, true);
  return tcc.p_count;
}

/// Entry point for the boost test executable
int main(int argc, char *argv[])
{
  using namespace precice;

  precice::syncMode = false;
  logging::setupLogging(); // first logging initalization, as early as possible
  utils::Parallel::initializeMPI(&argc, &argv);
  const auto rank = utils::Parallel::current()->rank();
  const auto size = utils::Parallel::current()->size();
  logging::setMPIRank(rank);

  if (size < 4) {
    if (rank == 0) {
      std::cerr << "ERROR: The tests require at least 4 MPI processes.\n";
    }
    return 2;
  }

  std::cout << "This test suite runs on rank " << rank << " of " << size << '\n';

  int       retCode  = boost::unit_test::unit_test_main(&init_unit_test, argc, argv);
  const int testsRan = countEnabledTests();

  // Override the return code if the slaves have nothing to test
  if ((testsRan == 0) && (rank != 0)) {
    retCode = EXIT_SUCCESS;
  }

  utils::MasterSlave::_communication = nullptr;
  utils::Parallel::finalizeMPI();
  return retCode;
}
