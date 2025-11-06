#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include "helpers.hpp"

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Serial)
BOOST_AUTO_TEST_SUITE(RBFParameterTuner)
PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(GaussianSupportRadiusPolynomial)
{
  /**
   * @brief Tests the Gaussian RBF mapping using the support radius option and a separate polynomial on linear data.
   *
   * The polynomial should explain the data almost completely leading to an error close to zero even before tuning.
   * Optimization will therefore stop before executing the bisection searching algorithm.
   */
  PRECICE_TEST();
  testRBFTuning(context.config(), context);
}

PRECICE_TEST_SETUP("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank))
BOOST_AUTO_TEST_CASE(GaussianSupportRadiusPUM)
{
  /**
   * @brief Tests the Gaussian RBF mapping using the support radius option on linear data. Optimization is only performed once.
   *
   * This test is separate from GaussianSupportRadiusPolynomial, to actually perform some bisection search iterations.
   */
  PRECICE_TEST();
  testRBFTuning(context.config(), context);
}


BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Serial
BOOST_AUTO_TEST_SUITE_END() // MappingRbfGaussian

#endif // PRECICE_NO_MPI
