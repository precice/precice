#include <Eigen/Core>
#include "../impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace cplscheme;

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

BOOST_AUTO_TEST_CASE(AbsoluteConvergenceMeasureTest)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  // Create convergence measure for Vector data
  double                                      convergenceLimit = 9.0;
  cplscheme::impl::AbsoluteConvergenceMeasure measure(convergenceLimit);

  // Create data sets for old state of data and new state of data
  Vector3d oldValues0(-2, -1, 0);
  Vector3d oldValues1(2, 3, 4);
  Vector3d oldValues2(3, 4, 5);
  Vector3d newValues(5, 6, 7);

  measure.measure(oldValues0, newValues);
  BOOST_TEST(not measure.isConvergence());

  measure.measure(oldValues1, newValues);
  BOOST_TEST(measure.isConvergence());

  measure.measure(oldValues2, newValues);
  BOOST_TEST(measure.isConvergence());
}

BOOST_AUTO_TEST_SUITE_END()
