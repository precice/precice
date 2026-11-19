#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(MappingRbfGaussian)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(GaussianSupportRadius)
{
  PRECICE_TEST();
  /**
   * @brief Tests the Gaussian rbf mapping using the support radius option
   *
   */
  testRBFMapping(context.config(), context);
}

PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(GaussianVectorialConsistent)
{
  PRECICE_TEST();
  /**
   * @brief Tests the Gaussian rbf mapping using vectorial data and a consistent mapping
   *
   */
  testRBFMappingVectorial(context.config(), context, false);
}

PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(GaussianVectorialConservative)
{
  PRECICE_TEST();
  /**
   * @brief Tests the Gaussian rbf mapping using vectorial data and a conservative mapping
   *
   */
  testRBFMappingVectorial(context.config(), context, true);
}
BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MappingRbfGaussian

#endif // PRECICE_NO_MPI
