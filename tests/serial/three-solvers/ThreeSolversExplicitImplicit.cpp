#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(ThreeSolvers)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank), "SolverThree"_on(1_rank))
BOOST_AUTO_TEST_CASE(ThreeSolversExplicitImplicit)
{
  PRECICE_TEST();
  std::string      config = context.config();
  std::vector<int> expectedCallsOfAdvance{30, 10, 30};
  runTestThreeSolvers(config, expectedCallsOfAdvance, context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // ThreeSolvers

#endif // PRECICE_NO_MPI
