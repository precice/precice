#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MapIfNecessary)
BOOST_AUTO_TEST_SUITE(TwoSolvers)
BOOST_AUTO_TEST_SUITE(WithoutSubsteps)

BOOST_AUTO_TEST_CASE(SerialImplicit)
{
  PRECICE_TEST("One"_on(1_rank), "Two"_on(1_rank));

  // 1: new end only (iterating)
  // 2: new end
  std::vector<int> readMappings{1, 1, 1, 1, 1, 0};
  // 1: we only map the sample from the end of the time window
  std::vector<int> writeMappings{1, 1, 1, 1, 1, 1};

  runTwoSolversMappingCountImplicit(context, readMappings, writeMappings);
}

BOOST_AUTO_TEST_SUITE_END() // WithoutSubsteps
BOOST_AUTO_TEST_SUITE_END() // TwoSolvers
BOOST_AUTO_TEST_SUITE_END() // MapIfNecessary
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
