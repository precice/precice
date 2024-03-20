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

BOOST_AUTO_TEST_CASE(SerialImplicit)
{
  PRECICE_TEST("A"_on(1_rank), "B"_on(1_rank), "C"_on(1_rank));
  // new data from A = new end
  // iterating B = new mid + new end
  // new tw B = new start + mid + end
  std::vector<int> readMappings{
      // initialize is not checked (should be start + mid + end from B and start + end from A (tw0))
      3, // iterating B + new data from A (tw1)
      2, // next tw: end from B
      3, // iterating B + final data from A (tw2)
      2, // next tw: end from B
      2, // iterating B
      0  // nothing from B (second)
  };
  // always maps mid + end per data
  std::vector<int> writeMappings{4, 4, 4, 4, 4, 4};

  runMultipleSolversMappingCount(context, readMappings, writeMappings);
}

BOOST_AUTO_TEST_SUITE_END() // MixedSubsteps
BOOST_AUTO_TEST_SUITE_END() // TwoSolvers
BOOST_AUTO_TEST_SUITE_END() // MapIfNecessary
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
