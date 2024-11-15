#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MapIfNecessary)
BOOST_AUTO_TEST_SUITE(ThreeSolvers)
BOOST_AUTO_TEST_SUITE(MixedSubsteps)

BOOST_AUTO_TEST_CASE(ParallelImplicit)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), "C"_on(1_rank));

  std::vector<int> readMappings{
      2, // iterating: mid and end from B
      2, // end from A + new start from B
      2, // iterating: mid and end from B
      2, // end from A + new start from B
      2, // iterating: mid and end from B
      2  // last end from A + last end from B
  };

  // When the coupling scheme A - B moves on to the next time window, C discards samples in the send data
  // namely MeshA:DataC.
  // While the itererative scheme is iterating, it still has this beginning timestamp and will be mapped in the writeMapping
  std::vector<int> writeMappings{4, 4, 4, 4, 4, 4};

  runMultipleSolversMappingCount(context, readMappings, writeMappings);
}

BOOST_AUTO_TEST_SUITE_END() // MixedSubsteps
BOOST_AUTO_TEST_SUITE_END() // TwoSolvers
BOOST_AUTO_TEST_SUITE_END() // MapIfNecessary
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
