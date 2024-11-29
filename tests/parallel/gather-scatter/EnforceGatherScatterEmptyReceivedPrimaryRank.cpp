#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_SUITE(GatherScatter)
PRECICE_TEST_SETUP("ParallelSolver"_on(2_ranks), "SerialSolver"_on(1_rank))
BOOST_AUTO_TEST_CASE(EnforceGatherScatterEmptyReceivedPrimaryRank)
{
  PRECICE_TEST();
  // Provided primary partition is not empty, but received primary partitionis empty
  runTestEnforceGatherScatter(std::vector<double>{0.0, 2.0, 0.0, 2.5}, context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel
BOOST_AUTO_TEST_SUITE_END() // GatherScatter

#endif // PRECICE_NO_MPI
