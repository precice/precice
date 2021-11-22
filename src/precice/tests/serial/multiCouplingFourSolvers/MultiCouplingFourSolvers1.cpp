#ifndef PRECICE_NO_MPI

#include <string>

#include "precice/tests/serial/multiCouplingFourSolvers/helpers.hpp"
#include "testing/Testing.hpp"
#include <fstream>

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_FIXTURE_TEST_SUITE(MultiCouplingFourSolvers, MultiCouplingFourSolversFixture)

BOOST_AUTO_TEST_CASE(MultiCouplingFourSolvers1)
{
  PRECICE_TEST("SolverA"_on(1_rank), "SolverB"_on(1_rank), "SolverC"_on(1_rank), "SolverD"_on(1_rank));
  const std::string configFile = _pathToTests + "multi-coupling-four-solver-1.xml";
  multiCouplingFourSolvers(configFile, context);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI