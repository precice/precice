#include "math/bspline.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

#include <Eigen/Core>

using namespace precice::math::bspline;
using namespace precice::testing;

BOOST_AUTO_TEST_SUITE(MathTests)
BOOST_AUTO_TEST_SUITE(BSpline)

BOOST_AUTO_TEST_CASE(TwoPointsLinear)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector2d ts;
  ts << 1, 2;
  Eigen::MatrixXd xs(3, 2);
  xs << 1, 2, 10, 20, 100, 200;

  // Limits
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 1.0), Eigen::Vector3d(1, 10, 100)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 2.0), Eigen::Vector3d(2, 20, 200)));

  // Midpoint
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 1.5), Eigen::Vector3d(1.5, 15, 150)));

  // Quarters
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 1.25), Eigen::Vector3d(1.25, 12.5, 125)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 1.75), Eigen::Vector3d(1.75, 17.5, 175)));
}

BOOST_AUTO_TEST_CASE(ThreePointsLinear)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d ts;
  ts << 0, 1, 2;
  Eigen::MatrixXd xs(3, 3);
  xs << 1, 2, 3, 10, 20, 30, 100, 200, 300;

  // Points
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 0.0), Eigen::Vector3d(1, 10, 100)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 1.0), Eigen::Vector3d(2, 20, 200)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 2.0), Eigen::Vector3d(3, 30, 300)));

  // Midpoints
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 0.5), Eigen::Vector3d(1.5, 15, 150)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 1.5), Eigen::Vector3d(2.5, 25, 250)));
}

BOOST_AUTO_TEST_CASE(ThreePointsLinearNonEquidistant)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d ts;
  ts << 0, 1, 3;
  Eigen::MatrixXd xs(3, 3);
  xs << 1, 2, 3, 10, 20, 30, 100, 200, 300;

  // Points
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 0.0), Eigen::Vector3d(1, 10, 100)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 1.0), Eigen::Vector3d(2, 20, 200)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 3.0), Eigen::Vector3d(3, 30, 300)));

  // Midpoints
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 0.5), Eigen::Vector3d(1.5, 15, 150)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 1, 2.0), Eigen::Vector3d(2.5, 25, 250), 1e-13));
}

BOOST_AUTO_TEST_CASE(ThreePointsQuadratic)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d ts;
  ts << 0, 1, 2;
  Eigen::MatrixXd xs(3, 3);
  xs << 1, 2, 3, 10, 20, 30, 100, 200, 300;

  // Points
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 0.0), Eigen::Vector3d(1, 10, 100), 1e-13));
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 1.0), Eigen::Vector3d(2, 20, 200)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 2.0), Eigen::Vector3d(3, 30, 300)));

  // Midpoints
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 0.5), Eigen::Vector3d(1.5, 15, 150)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 1.5), Eigen::Vector3d(2.5, 25, 250)));
}

BOOST_AUTO_TEST_CASE(ThreePointsQuadraticNonEquidistant)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d ts;
  ts << 0, 1, 3;
  Eigen::MatrixXd xs(3, 3);
  xs << 1, 2, 3, 10, 20, 30, 100, 200, 300;

  // Points
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 0.0), Eigen::Vector3d(1, 10, 100), 1e-13));
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 1.0), Eigen::Vector3d(2, 20, 200)));
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 3.0), Eigen::Vector3d(3, 30, 300)));

  // Midpoints
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 0.5), Eigen::Vector3d(1.54166, 15.41666, 154.16666), 1e-5));
  BOOST_TEST(equals(interpolateAt(ts, xs, 2, 2.0), Eigen::Vector3d(2.66666, 26.66666, 266.66666), 1e-5));
}

BOOST_AUTO_TEST_SUITE_END() // BSpline
BOOST_AUTO_TEST_SUITE_END() // Math
