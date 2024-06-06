#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MapIfNecessary)
BOOST_AUTO_TEST_SUITE(ThreeSolvers)
BOOST_AUTO_TEST_SUITE(WithSubsteps)

BOOST_AUTO_TEST_CASE(SerialImplicit)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), "C"_on(1_rank));

  std::vector<int> readMappings{
      // initialize not checked (tw0 of A)
      2, // iterating: mid + end from B
      3, // end from A (tw1) and mid + end from B
      2, // iterating: mid + end from B
      3, // final end from A (tw2) and mid + end from B
      2, // iterating: mid + end from B
      0  // nothing from B (second)
  };
  std::vector<int> writeMappings{4, 4, 4, 4, 4, 4};

  runMultipleSolversMappingCount(context, readMappings, writeMappings);
}

BOOST_AUTO_TEST_SUITE_END() // WithSubsteps
BOOST_AUTO_TEST_SUITE_END() // TwoSolvers
BOOST_AUTO_TEST_SUITE_END() // MapIfNecessary
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
