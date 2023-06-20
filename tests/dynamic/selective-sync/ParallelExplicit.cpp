#ifndef PRECICE_NO_MPI

#include "../../serial/explicit/helpers.hpp"
#include "testing/Testing.hpp"

#include <precice/precice.hpp>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Dynamic)
BOOST_AUTO_TEST_SUITE(SelectiveSync)
BOOST_AUTO_TEST_CASE(ParallelExplicit)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  runTestExplicit(context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Dynamic
BOOST_AUTO_TEST_SUITE_END() // SelectiveSync

#endif // PRECICE_NO_MPI
