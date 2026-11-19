#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/data/test_case.hpp>
#include <precice/precice.hpp>
#include "../helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(QuasiNewton)
BOOST_AUTO_TEST_SUITE(Serial)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_DATA_TEST_CASE(IMVJZeroUpdate, boost::unit_test::data::make({true, false}), includeSecondaryData)
{
  PRECICE_TEST();
  // serial coupling, IQN-IMVJ acceleration, to test the warning information to be given for zero update in QN-method
  runTestQN(includeSecondaryData, context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // QuasiNewton
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
