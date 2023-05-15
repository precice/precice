#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/Participant.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingRbfGaussian)
BOOST_AUTO_TEST_CASE(GaussianShapeParameter)
{
  /**
 * @brief Tests the Gaussian rbf mapping using the shape parameter option
 *
 */
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  testRBFMapping(context.config(), context);
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MappingRbfGaussian

#endif // PRECICE_NO_MPI
