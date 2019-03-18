#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_parameters.hpp>
#include <boost/filesystem.hpp>
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"
#include "utils/EventUtils.hpp"
#include "logging/LogConfiguration.hpp"
#include <iostream>

namespace precice {
extern bool testMode;
extern bool syncMode;
}


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
  
  auto & master_suite = framework::master_test_suite();
  master_suite.p_name.value = "preCICE Tests";

  auto logConfigs = logging::readLogConfFile("log.conf");
  
  if (logConfigs.empty()) { // nothing has been read from log.conf
    #if BOOST_VERSION == 106900
    std::cerr << "Boost 1.69 get log_level is broken, preCICE log level set to debug." << std::endl;
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
    
    logConfigs.push_back(config);
  }

  // Initialize either with empty logConfig -> default, configs that are read from file
  // or from the Boost Test log level.
  logging::setupLogging(logConfigs);
  logging::lockConf();
  

  // Sets the default tolerance for floating point comparisions
  // Can be overwritten on a per-test or per-suite basis using decators
  // boost::unit_test::decorator::collector::instance() * boost::unit_test::tolerance(0.001);
  * tolerance(1e-9); // Stores the decorator in the collector singleton
  #if BOOST_VERSION < 106900
  decorator::collector::instance().store_in(master_suite);
  #else
  decorator::collector_t::instance().store_in(master_suite);
  #endif
  
  return true;
}

/// Entry point for the boost test executable
int main(int argc, char* argv[])
{
  using namespace precice;

  precice::testMode = true;
  precice::syncMode = false;
  logging::setupLogging(); // first logging initalization, as early as possible
  utils::Parallel::initializeMPI(&argc, &argv);
  logging::setMPIRank(utils::Parallel::getProcessRank());
  utils::Petsc::initialize(&argc, &argv);
  utils::EventRegistry::instance().initialize("precice-Tests", "", utils::Parallel::getGlobalCommunicator());
    
  if (utils::Parallel::getCommunicatorSize() < 4) {
    if (utils::Parallel::getProcessRank() == 0)
      std::cerr << "Running tests on less than four processors. Not all tests are executed." << std::endl;
  }
  if (utils::Parallel::getCommunicatorSize() > 4) {
    if (utils::Parallel::getProcessRank() == 0)
      std::cerr << "Running tests on more than 4 processors is not supported. Aborting." << std::endl;
    std::exit(-1);
  }

  int retCode = boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

  utils::EventRegistry::instance().finalize();
  utils::Petsc::finalize();
  utils::Parallel::finalizeMPI();
  return retCode;
}
