#include "math/Bspline.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

#include <Eigen/Core>
#include <boost/test/data/test_case.hpp>

using namespace precice::testing;

BOOST_AUTO_TEST_SUITE(MathTests)
BOOST_AUTO_TEST_SUITE(BSpline)

void duplicateEndSample(Eigen::VectorXd &ts, Eigen::MatrixXd &xs, double epsilon = 1e-13)
{
  auto ntimestamps = ts.size();
  ts.conservativeResize(ntimestamps + 1);
  ts(ntimestamps) = ts(ntimestamps - 1);
  ts(ntimestamps - 1) -= epsilon;

  auto ndata = xs.rows();
  xs.conservativeResize(ndata, ntimestamps + 1);
  xs.col(ntimestamps) = xs.col(ntimestamps - 1);
}

BOOST_DATA_TEST_CASE(TwoPointsLinear,
                     boost::unit_test::data::make({false, true}),
                     duplicate)
{
  PRECICE_TEST(1_rank);
  Eigen::VectorXd ts = Eigen::Vector2d(1, 2);
  Eigen::MatrixXd xs(3, 2);
  xs << 1, 2, 10, 20, 100, 200;

  if (duplicate)
    duplicateEndSample(ts, xs);
  precice::math::Bspline bspline(ts, xs, 1);

  // Limits
  BOOST_TEST(equals(bspline.interpolateAt(1.0), Eigen::Vector3d(1, 10, 100)));
  BOOST_TEST(equals(bspline.interpolateAt(2.0), Eigen::Vector3d(2, 20, 200)));

  // Midpoint
  BOOST_TEST(equals(bspline.interpolateAt(1.5), Eigen::Vector3d(1.5, 15, 150)));

  // Quarters
  BOOST_TEST(equals(bspline.interpolateAt(1.25), Eigen::Vector3d(1.25, 12.5, 125)));
  BOOST_TEST(equals(bspline.interpolateAt(1.75), Eigen::Vector3d(1.75, 17.5, 175)));
}

BOOST_DATA_TEST_CASE(ThreePointsLinear,
                     boost::unit_test::data::make({false, true}),
                     duplicate)
{
  PRECICE_TEST(1_rank);
  Eigen::VectorXd ts = Eigen::Vector3d(0, 1, 2);
  Eigen::MatrixXd xs(3, 3);
  xs << 1, 2, 3, 10, 20, 30, 100, 200, 300;

  if (duplicate)
    duplicateEndSample(ts, xs);
  precice::math::Bspline bspline(ts, xs, 1);

  // Points
  BOOST_TEST(equals(bspline.interpolateAt(0.0), Eigen::Vector3d(1, 10, 100)));
  BOOST_TEST(equals(bspline.interpolateAt(1.0), Eigen::Vector3d(2, 20, 200)));
  BOOST_TEST(equals(bspline.interpolateAt(2.0), Eigen::Vector3d(3, 30, 300)));

  // Midpoints
  BOOST_TEST(equals(bspline.interpolateAt(0.5), Eigen::Vector3d(1.5, 15, 150)));
  BOOST_TEST(equals(bspline.interpolateAt(1.5), Eigen::Vector3d(2.5, 25, 250)));
}

BOOST_DATA_TEST_CASE(ThreePointsLinearNonEquidistant,
                     boost::unit_test::data::make({false, true}),
                     duplicate)
{
  Eigen::VectorXd ts = Eigen::Vector3d(0, 1, 3);
  Eigen::MatrixXd xs(3, 3);
  xs << 1, 2, 3, 10, 20, 30, 100, 200, 300;

  if (duplicate)
    duplicateEndSample(ts, xs);
  precice::math::Bspline bspline(ts, xs, 1);
  // Points
  BOOST_TEST(equals(bspline.interpolateAt(0.0), Eigen::Vector3d(1, 10, 100)));
  BOOST_TEST(equals(bspline.interpolateAt(1.0), Eigen::Vector3d(2, 20, 200)));
  BOOST_TEST(equals(bspline.interpolateAt(3.0), Eigen::Vector3d(3, 30, 300)));

  // Midpoints
  BOOST_TEST(equals(bspline.interpolateAt(0.5), Eigen::Vector3d(1.5, 15, 150)));
  BOOST_TEST(equals(bspline.interpolateAt(2.0), Eigen::Vector3d(2.5, 25, 250), 1e-13));
}

BOOST_DATA_TEST_CASE(ThreePointsQuadratic,
                     boost::unit_test::data::make({false, true}),
                     duplicate)
{
  PRECICE_TEST(1_rank);
  Eigen::VectorXd ts = Eigen::Vector3d(0, 1, 2);
  Eigen::MatrixXd xs(3, 3);
  xs << 1, 2, 3, 10, 20, 30, 100, 200, 300;

  if (duplicate)
    duplicateEndSample(ts, xs);
  precice::math::Bspline bspline(ts, xs, 2);
  // Points
  BOOST_TEST(equals(bspline.interpolateAt(0.0), Eigen::Vector3d(1, 10, 100), 1e-13));
  BOOST_TEST(equals(bspline.interpolateAt(1.0), Eigen::Vector3d(2, 20, 200)));
  BOOST_TEST(equals(bspline.interpolateAt(2.0), Eigen::Vector3d(3, 30, 300)));

  // Midpoints
  BOOST_TEST(equals(bspline.interpolateAt(0.5), Eigen::Vector3d(1.5, 15, 150)));
  BOOST_TEST(equals(bspline.interpolateAt(1.5), Eigen::Vector3d(2.5, 25, 250)));
}

BOOST_DATA_TEST_CASE(ThreePointsQuadraticNonEquidistant,
                     boost::unit_test::data::make({false, true}),
                     duplicate)
{
  PRECICE_TEST(1_rank);
  Eigen::VectorXd ts = Eigen::Vector3d(0, 1, 3);
  Eigen::MatrixXd xs(3, 3);
  xs << 1, 2, 3, 10, 20, 30, 100, 200, 300;

  if (duplicate)
    duplicateEndSample(ts, xs);
  precice::math::Bspline bspline(ts, xs, 2);
  // Points
  BOOST_TEST(equals(bspline.interpolateAt(0.0), Eigen::Vector3d(1, 10, 100), 1e-13));
  BOOST_TEST(equals(bspline.interpolateAt(1.0), Eigen::Vector3d(2, 20, 200)));
  BOOST_TEST(equals(bspline.interpolateAt(3.0), Eigen::Vector3d(3, 30, 300)));

  // Midpoints
  BOOST_TEST(equals(bspline.interpolateAt(0.5), Eigen::Vector3d(37.0 / 24, 370.0 / 24, 3700.0 / 24), 1e-13));
  BOOST_TEST(equals(bspline.interpolateAt(2.0), Eigen::Vector3d(8.0 / 3, 80.0 / 3, 800.0 / 3), 1e-13));
}

BOOST_AUTO_TEST_SUITE_END() // BSpline
BOOST_AUTO_TEST_SUITE_END() // Math
