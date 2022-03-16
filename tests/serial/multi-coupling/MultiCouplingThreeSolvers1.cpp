#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MultiCoupling)
BOOST_AUTO_TEST_CASE(MultiCouplingThreeSolvers1)
{
  PRECICE_TEST("SolverA"_on(1_rank), "SolverB"_on(1_rank), "SolverC"_on(1_rank));
  const std::string configFile = context.config();
  multiCouplingThreeSolvers(configFile, context);
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MultiCoupling

#endif // PRECICE_NO_MPI
