#ifndef PRECICE_NO_MPI

#include <fstream>
#include <string>
#include "helpers.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_CASE(MultiCouplingFourSolvers1)
{
  PRECICE_TEST("SolverA"_on(1_rank), "SolverB"_on(1_rank), "SolverC"_on(1_rank), "SolverD"_on(1_rank));
  const std::string configFile = context.config();
  multiCouplingFourSolvers(configFile, context);
}

BOOST_AUTO_TEST_SUITE_END() // PreciceTests
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
