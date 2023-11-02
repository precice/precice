#include <Eigen/Core>
#include "../impl/AbsoluteOrRelativeConvergenceMeasure.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace cplscheme;

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

BOOST_AUTO_TEST_CASE(AbsoluteOrRelativeConvergenceMeasureTest)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  // Create convergence measure for Vector data
  double                                                convergenceLimit        = 1.0;
  double                                                convergenceLimitPercent = 0.1;
  cplscheme::impl::AbsoluteOrRelativeConvergenceMeasure measure(convergenceLimit, convergenceLimitPercent);

  // Create data sets for old state of data and new state of data
  Vector3d oldValues0(1, 1, 1);
  Vector3d oldValues1(2, 2, 2);
  Vector3d oldValues2(2.6, 2.6, 2.6);
  Vector3d oldValues3(2.9, 2.9, 2.9);
  Vector3d newValues(3, 3, 3);

  measure.measure(oldValues0, newValues);
  BOOST_TEST(not measure.isConvergence());

  measure.measure(oldValues1, newValues);
  BOOST_TEST(not measure.isConvergence());

  measure.measure(oldValues2, newValues);
  BOOST_TEST(measure.isConvergence());

  measure.measure(oldValues3, newValues);
  BOOST_TEST(measure.isConvergence());
}

BOOST_AUTO_TEST_SUITE_END()
