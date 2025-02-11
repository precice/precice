#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include "../helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(QuasiNewton)
BOOST_AUTO_TEST_SUITE(Serial)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(TestQN11)
{
  PRECICE_TEST();
  // serial coupling,no subcycling but with interpolation, IQN-ILS
  runTestQNWithWaveforms(context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // QuasiNewton
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
