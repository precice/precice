#include "testing/Testing.hpp"
#include "utils/Parallel.hpp"

using namespace precice;

BOOST_AUTO_TEST_SUITE(TestingTests) // Use name of the module, e.g. subdirectory below src/, suffixed with Tests

BOOST_AUTO_TEST_SUITE(Examples) // If your file contains multiple tests, put them in a test suite

/// This tests runs indepentently on all processors
BOOST_AUTO_TEST_CASE(SingleProcessor)
{
  // Do not use DEBUG, TRACE, INFO calls inside tests
  
  BOOST_TEST(0 == 0); // Always use BOOST_TEST
}

/// Test with a modified numerical tolerance
BOOST_AUTO_TEST_CASE(NumericalTolerance,
                     * boost::unit_test::tolerance(1e-4))
{
  // Default tolerance is 1e-9, it can be changed for the entire case or even suite
  // using the decorator above
  BOOST_TEST(1.0 == 1.0001);

  // Or on a per test basis
  BOOST_TEST(1.0 == 1.01, boost::test_tools::tolerance(0.1));
}

/// Use testing::Deleted to unconditionally delete the test
BOOST_AUTO_TEST_CASE(Deleted,
                     * testing::Deleted())
{
  BOOST_TEST(false);
}


/// Tests that runs on 4 processors.
/*
 * If less than 4 procs are available, the test is deleted, if more are availablle, procs > 4 are deleted
 */
BOOST_AUTO_TEST_CASE(FourProcTests,
                     * testing::OnSize(4))
{
  BOOST_TEST(utils::Parallel::getCommunicatorSize() == 4);
}

/// Tests that runs on 2 processors.
/*
 * This case is trickier than with 4 procs, because we need to restrict the global communicator on all
 * ranks first, and then test if we execute at the correct ranks.
 */
BOOST_AUTO_TEST_CASE(TwoProcTests,
                     * testing::MinRanks(2)
                     * boost::unit_test::fixture<testing::MPICommRestrictFixture>(std::vector<int>({0, 1})))
{
  if (utils::Parallel::getCommunicatorSize() != 2)
    return;
  
  // Put your test code here
  
}


BOOST_AUTO_TEST_SUITE_END() // Examples
BOOST_AUTO_TEST_SUITE_END() // Testing
