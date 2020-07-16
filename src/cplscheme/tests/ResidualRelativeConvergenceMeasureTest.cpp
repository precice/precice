#include <Eigen/Core>
#include "../impl/ResidualRelativeConvergenceMeasure.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

BOOST_AUTO_TEST_CASE(ResidualRelativeConvergenceMeasureTest)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  double                                                       convergenceLimit = 0.1; // 10%
  precice::cplscheme::impl::ResidualRelativeConvergenceMeasure measure(convergenceLimit);

  // Create data sets for old state of data and new state of data
  Vector3d oldValues0(2.9, 2.9, 2.9);
  Vector3d oldValues1(2.95, 2.95, 2.95);
  Vector3d oldValues2(2.991, 2.991, 2.991);
  Vector3d newValues(3, 3, 3);

  // define initial residual, which we want to decrease in the following
  measure.measure(oldValues0, newValues);
  BOOST_TEST(not measure.isConvergence());

  measure.measure(oldValues1, newValues);
  BOOST_TEST(not measure.isConvergence());

  measure.measure(oldValues2, newValues);
  BOOST_TEST(measure.isConvergence());
}

BOOST_AUTO_TEST_SUITE_END()
