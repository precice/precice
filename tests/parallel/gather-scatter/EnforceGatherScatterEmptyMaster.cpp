#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(GatherScatter)
BOOST_AUTO_TEST_CASE(EnforceGatherScatterEmptyMaster)
{
  // Test case for an enforced gather scatter communication, where the partition
  // on the master rank is empty (recieved and provided). See issue #1013 for details.
  PRECICE_TEST("ParallelSolver"_on(2_ranks), "SerialSolver"_on(1_rank));
  // Provided master partition is empty and received master partition is empty
  runTestEnforceGatherScatter(std::vector<double>{}, context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // GatherScatter

#endif // PRECICE_NO_MPI
