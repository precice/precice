#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <boost/test/data/test_case.hpp>
#include <precice/precice.hpp>
#include "../helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(QuasiNewton)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_DATA_TEST_CASE(TestQN2, boost::unit_test::data::make({true, false}), includeSecondaryData)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  // parallel coupling, IQN-ILS, strict QR2 filter
  runTestQN(includeSecondaryData, context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // QuasiNewton
BOOST_AUTO_TEST_SUITE_END() // Serial

#endif // PRECICE_NO_MPI
