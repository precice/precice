#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MixedTimeWindowSizes)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_CASE(SerialParallel)
{
  PRECICE_TEST("Left"_on(1_rank), "Center"_on(1_rank), "Right"_on(1_rank));

  precice::testing::testThreeSolverStepOnlyExplicit(context);
}

BOOST_AUTO_TEST_SUITE_END() // Explicit
BOOST_AUTO_TEST_SUITE_END() // MixedTimeWindowSizes
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
