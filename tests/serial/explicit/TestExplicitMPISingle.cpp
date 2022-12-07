#ifndef PRECICE_NO_MPI

#include <precice/precice.hpp>
#include <vector>
#include "helpers.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(Explicit)
BOOST_AUTO_TEST_CASE(TestExplicitMPISingle)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  runTestExplicit(context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Explicit

#endif // PRECICE_NO_MPI
