#ifndef PRECICE_NO_MPI

#include "../helper.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MixedTimeWindowSizes)
BOOST_AUTO_TEST_SUITE(Implicit)
BOOST_AUTO_TEST_CASE(ParallelSerial)
{
  PRECICE_TEST("Left"_on(1_rank), "Center"_on(1_rank), "Right"_on(1_rank));

  precice::testing::testThreeSolverStepOnlyImplicit(context);
}

BOOST_AUTO_TEST_SUITE_END() // Implicit
BOOST_AUTO_TEST_SUITE_END() // MixedTimeWindowSizes
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration
#endif                      // PRECICE_NO_MPI
