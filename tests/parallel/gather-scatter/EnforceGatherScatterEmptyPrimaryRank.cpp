#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(GatherScatter)
PRECICE_TEST_SETUP("ParallelSolver"_on(2_ranks), "SerialSolver"_on(1_rank))
BOOST_AUTO_TEST_CASE(EnforceGatherScatterEmptyPrimaryRank)
{
  PRECICE_TEST();
  // Test case for an enforced gather scatter communication, where the partition
  // on the primary rank is empty (received and provided). See issue #1013 for details.
  // Provided primary partition is empty and received primary partition is empty
  runTestEnforceGatherScatter(std::vector<double>{}, context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // GatherScatter

#endif // PRECICE_NO_MPI
