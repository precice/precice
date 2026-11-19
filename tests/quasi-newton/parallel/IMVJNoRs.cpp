#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/data/test_case.hpp>
#include <precice/precice.hpp>
#include "../helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(QuasiNewton)
BOOST_AUTO_TEST_SUITE(Parallel)
PRECICE_TEST_SETUP("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks))
BOOST_DATA_TEST_CASE(IMVJNoRs, boost::unit_test::data::make({true, false}), includeSecondaryData)
{
  PRECICE_TEST();
  // serial coupling, IQN-IMVJ acceleration, to test `no-restart` method for restart;
  runTestQN(includeSecondaryData, context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // QuasiNewton
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI
