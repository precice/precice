#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include "../helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(QuasiNewton)
BOOST_AUTO_TEST_SUITE(Parallel)
PRECICE_TEST_SETUP("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks))
BOOST_AUTO_TEST_CASE(TestQN1EmptyPartition)
{
  PRECICE_TEST();
  // serial coupling, IQN-ILS, strict QR2 filter
  runTestQNEmptyPartition(context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // QuasiNewton
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI
