#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MapIfNecessary)
BOOST_AUTO_TEST_SUITE(ThreeSolvers)
BOOST_AUTO_TEST_SUITE(WithoutSubsteps)

PRECICE_TEST_SETUP("A"_on(1_rank), "B"_on(1_rank), "C"_on(1_rank))
BOOST_AUTO_TEST_CASE(SerialImplicit)
{
  PRECICE_TEST();

  std::vector<int> readMappings{
      // initialized not checked (tw0 of A)
      2, // iterate: mid + end of B
      3, // end of A (tw1) and new mid + end of B
      2, // iterate: mid + end of B
      3, // final end of A (tw2) and new mid + end of B
      2, // iterate: new mid + end of B
      0  // nothing from B (second)
  };
  std::vector<int> writeMappings{2, 2, 2, 2, 2, 2};

  runMultipleSolversMappingCount(context, readMappings, writeMappings);
}

BOOST_AUTO_TEST_SUITE_END() // WithoutSubsteps
BOOST_AUTO_TEST_SUITE_END() // TwoSolvers
BOOST_AUTO_TEST_SUITE_END() // MapIfNecessary
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
