// Setup boost test
// Disable the auto generation of main()
#define BOOST_TEST_NO_MAIN
// Specify the overall name of test framework
#define BOOST_TEST_MODULE "preCICE Tests"

#include <boost/test/tools/fpc_tolerance.hpp>
#include <boost/test/tree/test_case_counter.hpp>
#include <boost/test/tree/traverse.hpp>
#include <boost/test/unit_test.hpp>

#include <filesystem>
#include <iostream>
#include <string>

#include "com/SharedPointer.hpp"
#include "logging/LogConfiguration.hpp"
#include "utils/Ginkgo.hpp"
#include "utils/IntraComm.hpp"
#include "utils/Parallel.hpp"

namespace precice {
extern bool syncMode;
} // namespace precice

int countEnabledTests()
{
  using namespace boost::unit_test;
  test_case_counter tcc;
  traverse_test_tree(framework::master_test_suite(), tcc, true);
  return tcc.p_count;
}

class test_case_printer : public boost::unit_test::test_tree_visitor {
private:
  std::vector<std::string> prefix;

  bool test_suite_start(boost::unit_test::test_suite const &ts) override
  {
    if (ts.p_type_name == "suite") {
      prefix.push_back(ts.p_name);
    }
    return test_tree_visitor::visit((boost::unit_test::test_unit const &) ts);
  }

  void test_suite_finish(boost::unit_test::test_suite const &ts) override
  {
    if (ts.p_type_name == "suite") {
      prefix.pop_back();
    }
  }

  void visit(boost::unit_test::test_case const &tc) override
  {
    for (const auto &p : prefix) {
      std::cout << p << '/';
    }
    std::cout << tc.p_name << '\n';
  }
};

void printTestList()
{
  using namespace boost::unit_test;
  test_case_printer tcp;
  traverse_test_tree(framework::master_test_suite(), tcp, true);
}

void setupTolerance()
{
  static constexpr double tolerance = 1e-9;

  boost::test_tools::fpc_tolerance<double>() = tolerance;
  boost::test_tools::fpc_tolerance<float>()  = tolerance;
}

void removeStaleRunDirectory()
{
  namespace fs = std::filesystem;
  fs::path runDir("precice-run");
  if (fs::exists(runDir) && fs::is_directory(runDir) && !fs::is_empty(runDir)) {
    std::cout << "Removing a non-empty precice-run directory from a previously failing test.\n";
    fs::remove_all(runDir);
  }
}

/// Entry point for the boost test executable
int main(int argc, char *argv[])
{
  using namespace precice;

  precice::syncMode = false;
  utils::Parallel::initializeTestingMPI(&argc, &argv);
  const auto rank = utils::Parallel::current()->rank();
  const auto size = utils::Parallel::current()->size();
  logging::setMPIRank(rank);

  // Handle custom printing
  if (argc == 2 && std::string(argv[1]) == "--list_units") {
    if (rank == 0) {
      printTestList();
    }
    utils::Parallel::finalizeTestingMPI();
    return 0;
  }

  // Handle not enough MPI ranks
  if (size < 4 && argc < 2) {
    if (rank == 0) {
      std::cerr << "ERROR: The tests require at least 4 MPI processes. Please use \"mpirun -np 4 ./testprecice\" or \"ctest\" to run the full testsuite. \n";
    }
    utils::Parallel::finalizeTestingMPI();
    return 2;
  }

  std::cout << "This test suite runs on rank " << rank << " of " << size << '\n';

  if (rank == 0) {
    removeStaleRunDirectory();
  }

  setupTolerance();
  int       retCode  = boost::unit_test::unit_test_main(&init_unit_test, argc, argv);
  const int testsRan = countEnabledTests();

  // Override the return code if the secondary ranks have nothing to test
  if ((testsRan == 0) && (rank != 0)) {
    retCode = EXIT_SUCCESS;
  }
  // Required for Kokkos, which doesn't allow to initialize multiple times, i.e.,
  // finalize and initialize can really only be called once
#ifndef PRECICE_NO_GINKGO
  precice::utils::Ginkgo::finalize();
#endif
  utils::IntraComm::getCommunication() = nullptr;
  utils::Parallel::finalizeTestingMPI();
  return retCode;
}
