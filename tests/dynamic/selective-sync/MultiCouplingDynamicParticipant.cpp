#ifndef PRECICE_NO_MPI

#include "../../serial/multi-coupling/helpers.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Dynamic)
BOOST_AUTO_TEST_SUITE(SelectiveSync)
BOOST_AUTO_TEST_CASE(MultiCouplingDynamicParticipant)
{
  PRECICE_TEST("SolverA"_on(1_rank), "SolverB"_on(1_rank), "SolverC"_on(1_rank));
  const std::string configFile = context.config();
  multiCouplingThreeSolvers(configFile, context);
}

BOOST_AUTO_TEST_SUITE_END() // MultiCoupling
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // Integration

#endif // PRECICE_NO_MPI
