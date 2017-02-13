#include <boost/test/unit_test.hpp>
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"
#include "logging/LogConfiguration.hpp"

/// Boost test Initialization function:
bool init_unit_test()
{
  using namespace boost::unit_test;
  auto & master_suite = framework::master_test_suite();
  
  master_suite.p_name.value = "preCICE Tests";

  // Sets the default tolerance for floating point comparisions
  // Can be overwritten on a per-test or per-suite basis using decators
  // boost::unit_test::decorator::collector::instance() * boost::unit_test::tolerance(0.001);
  * tolerance(1e-9); // Stores the decorator in the collector singleton
  decorator::collector::instance().store_in(master_suite);
  
  return true;
}

/// Entry point for the boost test executable
int main(int argc, char* argv[])
{
  using namespace precice;

  logging::setupLogging();
  utils::Parallel::initializeMPI(&argc, &argv);
  logging::setMPIRank(utils::Parallel::getProcessRank());
  utils::Petsc::initialize(&argc, &argv);
  
  int retCode = boost::unit_test::unit_test_main( &init_unit_test, argc, argv );

  utils::Petsc::finalize();
  utils::Parallel::finalizeMPI();
  return retCode;
}
