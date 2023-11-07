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
  double                                                absLimit1 = 1.0;
  double                                                relLimit1 = 0.05;
  double                                                absLimit2 = 0.2;
  double                                                relLimit2 = 0.2;
  cplscheme::impl::AbsoluteOrRelativeConvergenceMeasure measure1(absLimit1, relLimit1);
  cplscheme::impl::AbsoluteOrRelativeConvergenceMeasure measure2(absLimit2, relLimit2);

  // Create data sets for old state of data and new state of data
  Vector3d oldValues0(1, 1, 1);
  Vector3d oldValues1(2.6, 2.6, 2.6);
  Vector3d oldValues2(2.9, 2.9, 2.9);
  Vector3d newValues(3, 3, 3);

  measure1.measure(oldValues0, newValues);
  BOOST_TEST(not measure1.isConvergence());

  measure1.measure(oldValues1, newValues);
  BOOST_TEST(measure1.isConvergence());

  measure1.measure(oldValues2, newValues);
  BOOST_TEST(measure1.isConvergence());

  measure2.measure(oldValues0, newValues);
  BOOST_TEST(not measure2.isConvergence());

  measure2.measure(oldValues1, newValues);
  BOOST_TEST(measure2.isConvergence());

  measure2.measure(oldValues2, newValues);
  BOOST_TEST(measure2.isConvergence());
}

BOOST_AUTO_TEST_SUITE_END()
