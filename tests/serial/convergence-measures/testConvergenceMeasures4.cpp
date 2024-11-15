#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(ConvergenceMeasures)
BOOST_AUTO_TEST_CASE(testConvergenceMeasures4)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  std::vector<int> expectedIterations = {2, 2};
  testConvergenceMeasures(context.config(), context, expectedIterations);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // ConvergenceMeasures

#endif // PRECICE_NO_MPI
