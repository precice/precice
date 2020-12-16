#include <Eigen/Core>
#include "../geometry.hpp"
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "math/differences.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "utils/algorithm.hpp"

using namespace precice;
using namespace precice::math;

BOOST_AUTO_TEST_SUITE(MathTests)
BOOST_AUTO_TEST_SUITE(Geometry)

BOOST_AUTO_TEST_CASE(Collinear)
{
  PRECICE_TEST(1_rank);
  // 2D test setup
  Eigen::Vector2d a2D(0, 0);
  Eigen::Vector2d b2D(1, 1);
  Eigen::Vector2d collinearPoint2D(0.5, 0.5);
  Eigen::Vector2d notCollinearPoint2D(0.5, 0.6);

  // 3D test setup
  Eigen::Vector3d a3D(0, 0, 0);
  Eigen::Vector3d b3D(1, 1, 1);
  Eigen::Vector3d collinearPoint3D(0.5, 0.5, 0.5);
  Eigen::Vector3d notCollinearPoint3D(0.5, 0.6, 0.5);

  // 2D test validations
  BOOST_CHECK(geometry::collinear(a2D, b2D, collinearPoint2D));
  BOOST_CHECK(!geometry::collinear(a2D, b2D, notCollinearPoint2D));

  // 3D test validations
  BOOST_CHECK(geometry::collinear(a3D, b3D, collinearPoint3D));
  BOOST_CHECK(!geometry::collinear(a3D, b3D, notCollinearPoint3D));
}

BOOST_AUTO_TEST_CASE(TetraVolume,
                     *boost::unit_test::tolerance(1e-3))
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d a(1, 2, 3);
  Eigen::Vector3d b(3, 2, 1);
  Eigen::Vector3d c(4, 5, 6);
  Eigen::Vector3d d(6, 5, 4);
  BOOST_TEST(geometry::tetraVolume(a, b, c, d) == 0);

  a << 5, 0, 0;
  b << 0, -3, -3;
  c << 0, 3, 4;
  d << -1, -2, 6;
  BOOST_TEST(geometry::tetraVolume(a, b, c, d) == 38.6666);

  d << -1.47, -4.1, 8.3;
  BOOST_TEST(geometry::tetraVolume(a, b, c, d) == 62.1816);
}

BOOST_AUTO_TEST_CASE(Between)
{
  PRECICE_TEST(1_rank);
  for (int dim = 2; dim <= 3; dim++) {
    Eigen::VectorXd a(dim);
    Eigen::VectorXd b(dim);
    Eigen::VectorXd betweenPoint(dim);
    Eigen::VectorXd betweenLimitPoint(dim);
    Eigen::VectorXd collinearOutsidePoint(dim);
    Eigen::VectorXd outsidePoint(dim);
    if (dim == 2) {
      a << 0.0, 0.0;
      b << 1.0, 1.0;
      betweenPoint << 0.5, 0.5;
      betweenLimitPoint << 1.0, 1.0;
      collinearOutsidePoint << 2.0, 2.0;
      outsidePoint << 0.5, 0.4;
    } else {
      a << 0.0, 0.0, 0.0;
      b << 1.0, 1.0, 1.0;
      betweenPoint << 0.5, 0.5, 0.5;
      betweenLimitPoint << 1.0, 1.0, 1.0;
      collinearOutsidePoint << 2.0, 2.0, 2.0;
      outsidePoint << 0.5, 0.4, 0.5;
    }
    BOOST_CHECK(geometry::between(a, b, betweenPoint));
    BOOST_CHECK(geometry::between(a, b, betweenLimitPoint));
    BOOST_CHECK(!geometry::between(a, b, collinearOutsidePoint));
    BOOST_CHECK(!geometry::between(a, b, outsidePoint));
  }
}

BOOST_AUTO_TEST_CASE(TriangleArea)
{
  PRECICE_TEST(1_rank);
  { // 2D
    Eigen::Vector2d a, b, c;
    double          area;
    a << 0.0, 0.0;
    b << 1.0, 0.0;
    c << 0.0, 1.0;
    area = geometry::triangleArea(a, b, c);
    BOOST_TEST(area == 0.5);

    b << 0.0, 1.0;
    c << 1.0, 0.0;
    area = geometry::triangleArea(a, b, c);
    BOOST_TEST(area == -0.5);
  }
  { // 3D
    Eigen::Vector3d a, b, c;
    double          area;
    a << 0.0, 0.0, 0.0;
    b << 1.0, 0.0, 1.0;
    c << 1.0, 1.0, 1.0;
    area = geometry::triangleArea(a, b, c);
    BOOST_CHECK(area == 0.5 * sqrt(2.0));
  }
}

BOOST_AUTO_TEST_CASE(SegmentsIntersect)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector2d a(0, 0), b(1, 0), c(0.5, 0), d(0, 0.5);
  BOOST_CHECK(geometry::segmentsIntersect(a, b, c, d, true));
  BOOST_CHECK(!geometry::segmentsIntersect(a, b, c, d, false));

  c(1) = -0.2;
  BOOST_CHECK(geometry::segmentsIntersect(a, b, c, d, false));

  // Test case motivated from bug
  a << 0.23104429, 1.87753905;
  b << 0.22608634, 1.88114120;
  c << 0.23058985, 1.87733882;
  d << 0.23058985, 1.87846349;
  BOOST_CHECK(geometry::segmentsIntersect(a, b, c, d, false));

  // Another bug motivated test, T-intersection slightly above numerical accuracy
  a << 1.0, 1.0;
  b << 0.0, 1.0;
  c << 1.0 - (10.0 * NUMERICAL_ZERO_DIFFERENCE), 1.1;
  d << 1.0 - (10.0 * NUMERICAL_ZERO_DIFFERENCE), 0.9;
  BOOST_CHECK(geometry::segmentsIntersect(a, b, c, d, false));
  BOOST_CHECK(geometry::segmentsIntersect(c, d, a, b, false));
}

BOOST_AUTO_TEST_CASE(SegmentPlaneIntersection)
{
  PRECICE_TEST(1_rank);
  using Eigen::Vector3d;
  Vector3d planeNormal        = Vector3d::Constant(1.0);
  Vector3d pointOnPlane       = Vector3d::Constant(0.0);
  Vector3d firstSegmentPoint  = Vector3d::Constant(1.0);
  Vector3d secondSegmentPoint = Vector3d::Constant(-1.0);
  Vector3d intersectionPoint  = Vector3d::Constant(1.0);
  Vector3d expected           = Vector3d::Constant(0.0);

  // True intersection
  int result = geometry::segmentPlaneIntersection(
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint);
  BOOST_TEST(result == geometry::INTERSECTION);
  BOOST_CHECK(equals(intersectionPoint, expected));

  // Touching second segment vertex
  secondSegmentPoint = Vector3d::Constant(0.0);
  result             = geometry::segmentPlaneIntersection(
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint);
  BOOST_TEST(result == geometry::TOUCHING);
  BOOST_CHECK(equals(intersectionPoint, expected));

  // Touching first segment vertex
  firstSegmentPoint  = Vector3d::Constant(0.0);
  secondSegmentPoint = Vector3d::Constant(-1.0);
  result             = geometry::segmentPlaneIntersection(
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint);
  BOOST_TEST(result == geometry::TOUCHING);
  BOOST_CHECK(equals(intersectionPoint, expected));

  // Parallel segment with distance to plain
  firstSegmentPoint << 0.0, 0.0, -3.0;
  intersectionPoint << 1.0, 2.0, 3.0; // should not be modified
  expected = intersectionPoint;
  result   = geometry::segmentPlaneIntersection(
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint);
  BOOST_TEST(result == geometry::NO_INTERSECTION);
  BOOST_CHECK(equals(intersectionPoint, expected));

  // Parallel segment contained in plane
  firstSegmentPoint << 0.0, 0.0, 0.0;
  secondSegmentPoint << 1.0, 1.0, -2.0;
  result = geometry::segmentPlaneIntersection(
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint);
  BOOST_TEST(result == geometry::CONTAINED);
  BOOST_CHECK(equals(intersectionPoint, expected));

  // Segment ending before intersection
  firstSegmentPoint << -2.0, -2.0, -2.0;
  secondSegmentPoint << -1.0, -1.0, -1.0;
  result = geometry::segmentPlaneIntersection(
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint);
  BOOST_TEST(result == geometry::NO_INTERSECTION);
  BOOST_CHECK(equals(intersectionPoint, expected));

  // Segment ending before intersection (inversed segment points)
  firstSegmentPoint << -1.0, -1.0, -1.0;
  secondSegmentPoint << -2.0, -2.0, -2.0;
  result = geometry::segmentPlaneIntersection(
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint);
  BOOST_TEST(result == geometry::NO_INTERSECTION);
  BOOST_CHECK(equals(intersectionPoint, expected));
}

BOOST_AUTO_TEST_CASE(ProjectVector)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector3d vector3D(1.0, 2.0, 3.0);
  Eigen::Vector2d vector2D;
  Eigen::Vector2d vectorExpected(1.0, 2.0);

  vector2D = geometry::projectVector(vector3D, 2);
  BOOST_CHECK(equals(vector2D, vectorExpected));

  vector2D = geometry::projectVector(vector3D, 1);
  vectorExpected << 1.0, 3.0;
  BOOST_CHECK(equals(vector2D, vectorExpected));

  vector2D = geometry::projectVector(vector3D, 0);
  vectorExpected << 2.0, 3.0;
  BOOST_CHECK(equals(vector2D, vectorExpected));
}

BOOST_AUTO_TEST_CASE(ContainedInTriangle)
{
  PRECICE_TEST(1_rank);
  Eigen::Vector2d triangleVertex0(0.0, 0.0);
  Eigen::Vector2d triangleVertex1(1.0, 0.0);
  Eigen::Vector2d triangleVertex2(0.0, 1.0);
  Eigen::Vector2d point(0.25, 0.25);

  // Contained point
  int result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::CONTAINED);

  // Not contained points
  point << -1.0, -1.0;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::NOT_CONTAINED);
  point << 1.0, 1.0;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::NOT_CONTAINED);
  point << 2.0, 0.0;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::NOT_CONTAINED);
  point << 0.0, 2.0;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  // Touching points
  point << 0.0, 0.0;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::TOUCHING);
  point << 1.0, 0.0;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::TOUCHING);
  point << 0.0, 1.0;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::TOUCHING);
  point << 0.5, 0.0;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::TOUCHING);
  point << 0.0, 0.5;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::TOUCHING);
  point << 0.5, 0.5;
  result = geometry::containedInTriangle(
      triangleVertex0, triangleVertex1, triangleVertex2, point);
  BOOST_TEST(result == geometry::TOUCHING);
}

BOOST_AUTO_TEST_CASE(ContainedInHyperrectangle)
{
  PRECICE_TEST(1_rank);
  // 2D
  Eigen::Vector2d center2D(0, 0);
  Eigen::Vector2d sidelengths2D(1, 1);

  // Not contained 2D
  Eigen::Vector2d testPoint2D(2, 2);
  int             result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint2D << -2.0, -2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint2D << 2.0, -2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint2D << -2.0, 2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint2D << 2.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint2D << -2.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint2D << 0.0, 2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint2D << 0.0, -2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  // Contained 2D
  testPoint2D << 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint2D << 0.25, 0.25;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint2D << -0.25, 0.25;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint2D << 0.25, -0.25;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint2D << -0.25, -0.25;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint2D << 0.49999999999, 0.49999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint2D << -0.49999999999, -0.49999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::CONTAINED);

  // Touching 2D
  testPoint2D << 0.5, 0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << -0.5, 0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.5, -0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << -0.5, -0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.4999999999999999, 0.4999999999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.500000000000001, 0.500000000000001;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << -0.500000000000001, -0.500000000000001;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.5, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << -0.5, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.0, 0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.0, -0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.5, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << -0.5, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.0, 0.499999999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.0, -0.499999999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.0, 0.500000000000001;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.0, -0.500000000000001;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.499999999999999, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << -0.499999999999999, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << 0.500000000000001, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint2D << -0.500000000000001, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths2D, center2D, testPoint2D);
  BOOST_TEST(result == geometry::TOUCHING);

  // 3D
  Eigen::Vector3d center3D      = Eigen::Vector3d::Zero();
  Eigen::Vector3d sidelengths3D = Eigen::Vector3d::Constant(1.0);

  // Not contained 3D
  Eigen::Vector3d testPoint3D(2, 2, 2);
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint3D << -2.0, -2.0, -2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint3D << 2.0, -2.0, 2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint3D << -2.0, 2.0, 2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint3D << 2.0, 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint3D << -2.0, 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint3D << 0.0, 2.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint3D << 0.0, 0.0, 2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  testPoint3D << 0.0, 0.0, -2.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::NOT_CONTAINED);

  // Contained 3D
  testPoint3D << 0.0, 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint3D << 0.25, 0.25, 0.25;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint3D << -0.25, 0.25, 0.25;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint3D << 0.25, -0.25, 0.25;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint3D << -0.25, -0.25, 0.25;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint3D << 0.25, 0.25, -0.25;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint3D << 0.49999999999, 0.49999999999, 0.49999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::CONTAINED);

  testPoint3D << -0.49999999999, -0.49999999999, -0.49999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::CONTAINED);

  // Touching 3D
  testPoint3D << 0.5, 0.5, 0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << -0.5, 0.5, 0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.5, -0.5, 0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << -0.5, -0.5, -0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.4999999999999999, 0.4999999999999999, 0.4999999999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.500000000000001, 0.500000000000001, 0.500000000000001;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << -0.500000000000001, -0.500000000000001, -0.500000000000001;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.5, 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << -0.5, 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, 0.5, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, -0.5, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.5, 0.0, -0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << -0.5, 0.0, 0.5;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, 0.499999999999999, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, -0.499999999999999, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, 0.500000000000001, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, -0.500000000000001, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.499999999999999, 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << -0.499999999999999, 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.500000000000001, 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << -0.500000000000001, 0.0, 0.0;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, 0.0, 0.499999999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, 0.0, -0.499999999999999;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, 0.0, 0.500000000000001;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);

  testPoint3D << 0.0, 0.0, -0.500000000000001;
  result = geometry::containedInHyperrectangle(
      sidelengths3D, center3D, testPoint3D);
  BOOST_TEST(result == geometry::TOUCHING);
}

BOOST_AUTO_TEST_SUITE(Convexity)

BOOST_AUTO_TEST_CASE(ComputeUnitQuadConvexity)
{
  PRECICE_TEST(1_rank);
  int             dim = 3;
  Eigen::VectorXd coords0(dim);
  Eigen::VectorXd coords1(dim);
  Eigen::VectorXd coords2(dim);
  Eigen::VectorXd coords3(dim);
  coords0 << 0.0, 0.0, 0.0;
  coords1 << 1.0, 0.0, 0.0;
  coords2 << 1.0, 1.0, 0.0;
  coords3 << 0.0, 1.0, 0.0;

  auto vertexList = utils::make_array(coords0, coords1, coords2, coords3);
  auto result     = geometry::isConvexQuad(vertexList);

  BOOST_TEST(result.convex);
  BOOST_TEST(utils::unique_elements(result.vertexOrder));
  BOOST_TEST_MESSAGE("Vertex Order" << result.vertexOrder);
  BOOST_TEST(result.vertexOrder.at(0) == 0);
  BOOST_TEST(result.vertexOrder.at(1) == 3);
  BOOST_TEST(result.vertexOrder.at(2) == 2);
  BOOST_TEST(result.vertexOrder.at(3) == 1);
}

BOOST_AUTO_TEST_CASE(ComputeReversedUnitQuadConvexity)
{
  PRECICE_TEST(1_rank);
  int             dim = 3;
  Eigen::VectorXd coords0(dim);
  Eigen::VectorXd coords1(dim);
  Eigen::VectorXd coords2(dim);
  Eigen::VectorXd coords3(dim);
  coords0 << 0.0, 1.0, 0.0;
  coords1 << 1.0, 1.0, 0.0;
  coords2 << 1.0, 0.0, 0.0;
  coords3 << 0.0, 0.0, 0.0;

  auto vertexList = utils::make_array(coords0, coords1, coords2, coords3);
  auto result     = geometry::isConvexQuad(vertexList);

  BOOST_TEST(result.convex);
  BOOST_TEST(utils::unique_elements(result.vertexOrder));
  BOOST_TEST_MESSAGE("Vertex Order" << result.vertexOrder);
  BOOST_TEST(result.vertexOrder.at(0) == 0);
  BOOST_TEST(result.vertexOrder.at(1) == 3);
  BOOST_TEST(result.vertexOrder.at(2) == 2);
  BOOST_TEST(result.vertexOrder.at(3) == 1);
}

BOOST_AUTO_TEST_CASE(ComputeValidQuadConvexity)
{
  PRECICE_TEST(1_rank);
  int             dim = 3;
  Eigen::VectorXd coords0(dim);
  Eigen::VectorXd coords1(dim);
  Eigen::VectorXd coords2(dim);
  Eigen::VectorXd coords3(dim);
  coords0 << 0.5, 0.34, 0.0;
  coords1 << 0.62, 0.32, 0.0;
  coords2 << 0.6, 0.24, 0.0;
  coords3 << 0.3, 0.22, 0.0;

  auto vertexList = utils::make_array(coords0, coords1, coords2, coords3);
  auto result     = geometry::isConvexQuad(vertexList);

  BOOST_TEST(result.convex);
  BOOST_TEST(utils::unique_elements(result.vertexOrder));
  BOOST_TEST_MESSAGE("Vertex Order" << result.vertexOrder);
  BOOST_TEST(result.vertexOrder.at(0) == 3);
  BOOST_TEST(result.vertexOrder.at(1) == 2);
  BOOST_TEST(result.vertexOrder.at(2) == 1);
  BOOST_TEST(result.vertexOrder.at(3) == 0);
}

BOOST_AUTO_TEST_CASE(ComputeValidQuadConvexityWithOffPlane)
{
  PRECICE_TEST(1_rank);
  int             dim = 3;
  Eigen::VectorXd coords0(dim);
  Eigen::VectorXd coords1(dim);
  Eigen::VectorXd coords2(dim);
  Eigen::VectorXd coords3(dim);
  coords0 << 0.5, 0.34, 0.0;
  coords1 << 0.62, 0.32, 0.0;
  coords2 << 0.6, 0.24, 0.0;
  coords3 << 0.3, 0.22, 0.5;

  auto vertexList = utils::make_array(coords0, coords1, coords2, coords3);
  // This test should result in an error
  // auto result = geometry::isConvexQuad(vertexList);
  //
  // BOOST_TEST(result.convex);
  // BOOST_TEST(utils::unique_elements(result.vertexOrder));
  // BOOST_TEST_MESSAGE("Vertex Order" << result.vertexOrder);
  // BOOST_TEST(result.vertexOrder.at(0) == 3);
  // BOOST_TEST(result.vertexOrder.at(1) == 2);
  // BOOST_TEST(result.vertexOrder.at(2) == 1);
  // BOOST_TEST(result.vertexOrder.at(3) == 0);
}

BOOST_AUTO_TEST_CASE(ComputeInvalidUnitQuadConvexity)
{
  PRECICE_TEST(1_rank);
  int             dim = 3;
  Eigen::VectorXd coords0(dim);
  Eigen::VectorXd coords1(dim);
  Eigen::VectorXd coords2(dim);
  Eigen::VectorXd coords3(dim);
  coords2 << 1.0, 0.0, 0.0;
  coords1 << 1.0, 1.0, 0.0;
  coords0 << 0.5, 0.9, 0.0;
  coords3 << 0.0, 1.0, 0.0;

  auto vertexList = utils::make_array(coords0, coords1, coords2, coords3);
  auto result     = geometry::isConvexQuad(vertexList);

  BOOST_TEST(!result.convex);
}

BOOST_AUTO_TEST_CASE(ComputeInvalidQuadConvexity)
{
  PRECICE_TEST(1_rank);
  int             dim = 3;
  Eigen::VectorXd coords0(dim);
  Eigen::VectorXd coords1(dim);
  Eigen::VectorXd coords2(dim);
  Eigen::VectorXd coords3(dim);
  coords0 << 0.5, 0.34, 0.0;
  coords1 << 0.62, 0.32, 0.0;
  coords2 << 0.52, 0.31, 0.0;
  coords3 << 0.51, 0.22, 0.0;

  auto vertexList = utils::make_array(coords0, coords1, coords2, coords3);
  auto result     = geometry::isConvexQuad(vertexList);

  BOOST_TEST(!result.convex);
}

BOOST_AUTO_TEST_SUITE_END() // convexity

BOOST_AUTO_TEST_SUITE_END() // geometry

BOOST_AUTO_TEST_SUITE_END() // Math
