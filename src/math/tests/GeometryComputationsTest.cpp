#include "GeometryComputationsTest.hpp"
#include "../GeometryComputations.hpp"
#include "../Parallel.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::GeometryComputationsTest)

namespace precice {
namespace utils {
namespace tests {

logging::Logger GeometryComputationsTest::_log ("precice::math::GeometryComputationsTest");

GeometryComputationsTest:: GeometryComputationsTest ()
:
  TestCase ("math::GeometryComputationsTest")
{}

void GeometryComputationsTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod(testTriangleArea);
    testMethod(testCollinear);
    testMethod(testBetween);
    testMethod(testSegmentsIntersect);
    testMethod(testSegmentPlaneIntersection);
    testMethod(testProjectVector);
    testMethod(testContainedInTriangle);
    testMethod(testContainedInHyperrectangle);
  }
}

void GeometryComputationsTest:: testCollinear ()
{
  TRACE();
  // 2D test setup
  Eigen::Vector2d a2D (0, 0);
  Eigen::Vector2d b2D (1, 1);
  Eigen::Vector2d collinearPoint2D (0.5, 0.5);
  Eigen::Vector2d notCollinearPoint2D (0.5, 0.6);

  // 3D test setup
  Eigen::Vector3d a3D (0, 0, 0);
  Eigen::Vector3d b3D (1, 1, 1);
  Eigen::Vector3d collinearPoint3D (0.5, 0.5, 0.5);
  Eigen::Vector3d notCollinearPoint3D (0.5, 0.6, 0.5);

  // 2D test validations
  validate ( GeometryComputations::collinear(a2D, b2D, collinearPoint2D) );
  validate ( ! GeometryComputations::collinear(a2D, b2D, notCollinearPoint2D) );

  // 3D test validations
  validate ( GeometryComputations::collinear(a3D, b3D, collinearPoint3D) );
  validate ( ! GeometryComputations::collinear(a3D, b3D, notCollinearPoint3D) );
}

void GeometryComputationsTest::testTetraVolume()
{
  TRACE();
  Eigen::Vector3d a(1,2,3);
  Eigen::Vector3d b(3,2,1);
  Eigen::Vector3d c(4,5,6);
  Eigen::Vector3d d(6,5,4);
  validateEquals( GeometryComputations::tetraVolume(a, b, c, d), 0 );

  a << 5, 0, 0;
  b << 0, -3, -3;
  c << 0, 3, 4;
  d << -1, -2, 6;
  validateEquals( GeometryComputations::tetraVolume(a, b, c, d), 232 );

  d << -1.47, -4.1, 8.3;
  validateEquals( GeometryComputations::tetraVolume(a, b, c, d), 373.3 );
}

void GeometryComputationsTest:: testBetween ()
{
  TRACE();
  for ( int dim=2; dim <= 3; dim++ ){
    Eigen::VectorXd a(dim);
    Eigen::VectorXd b(dim);
    Eigen::VectorXd betweenPoint(dim);
    Eigen::VectorXd betweenLimitPoint(dim);
    Eigen::VectorXd collinearOutsidePoint(dim);
    Eigen::VectorXd outsidePoint(dim);
    if ( dim == 2 ){
      a << 0.0, 0.0;
      b << 1.0, 1.0;
      betweenPoint << 0.5, 0.5;
      betweenLimitPoint << 1.0, 1.0;
      collinearOutsidePoint << 2.0, 2.0;
      outsidePoint << 0.5, 0.4;
    }
    else {
      a << 0.0, 0.0, 0.0;
      b << 1.0, 1.0, 1.0;
      betweenPoint << 0.5, 0.5, 0.5;
      betweenLimitPoint << 1.0, 1.0, 1.0;
      collinearOutsidePoint << 2.0, 2.0, 2.0;
      outsidePoint << 0.5, 0.4, 0.5;
    }
    validate(math::GeometryComputations::between(a, b, betweenPoint));
    validate(math::GeometryComputations::between(a, b, betweenLimitPoint));
    validate(! math::GeometryComputations::between(a, b, collinearOutsidePoint));
    validate(! math::GeometryComputations::between(a, b, outsidePoint));
  }
}

void GeometryComputationsTest:: testTriangleArea ()
{
  TRACE()
  { // 2D
    Eigen::Vector2d a, b, c;
    double area;
    a << 0.0, 0.0;
    b << 1.0, 0.0;
    c << 0.0, 1.0;
    area = math::GeometryComputations::triangleArea(a, b, c);
    validateNumericalEquals (area, 0.5);

    b << 0.0, 1.0;
    c << 1.0, 0.0;
    area = math::GeometryComputations::triangleArea(a, b, c);
    validateNumericalEquals (area, -0.5);
  }
  { // 3D
    Eigen::Vector3d a, b, c;
    double area;
    a << 0.0, 0.0, 0.0;
    b << 1.0, 0.0, 1.0;
    c << 1.0, 1.0, 1.0;
    area = math::GeometryComputations::triangleArea(a, b, c);
    validateNumericalEquals (area, 0.5 * sqrt(2.0));
  }
}

void GeometryComputationsTest:: testSegmentsIntersect ()
{
  TRACE();
  typedef GeometryComputations GeoComp;
  Eigen::Vector2d a(0,0), b(1,0), c(0.5,0), d(0,0.5);
  validate ( GeoComp::segmentsIntersect(a,b,c,d,true));
  validate ( ! GeoComp::segmentsIntersect(a,b,c,d,false));

  c(1) = -0.2;
  validate ( GeoComp::segmentsIntersect(a,b,c,d,false));

  // Test case motivated from bug
  a << 0.23104429, 1.87753905;
  b << 0.22608634, 1.88114120;
  c << 0.23058985, 1.87733882;
  d << 0.23058985, 1.87846349;
  validate ( GeoComp::segmentsIntersect(a,b,c,d,false));

  // Another bug motivated test, T-intersection slightly above numerical accuracy
  a << 1.0, 1.0;
  b << 0.0, 1.0;
  c << 1.0 - (10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE), 1.1;
  d << 1.0 - (10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE), 0.9;
  validate ( GeoComp::segmentsIntersect(a,b,c,d,false) );
  validate ( GeoComp::segmentsIntersect(c,d,a,b,false) );
}

void GeometryComputationsTest:: testSegmentPlaneIntersection ()
{
  TRACE();
  using Eigen::Vector3d;
  Vector3d planeNormal = Vector3d::Constant ( 1.0 );
  Vector3d pointOnPlane = Vector3d::Constant ( 0.0 );
  Vector3d firstSegmentPoint = Vector3d::Constant ( 1.0 );
  Vector3d secondSegmentPoint = Vector3d::Constant ( -1.0 );
  Vector3d intersectionPoint = Vector3d::Constant ( 1.0 );
  Vector3d expected = Vector3d::Constant ( 0.0 );

  // True intersection
  int result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::INTERSECTION );
  validate ( math::equals(intersectionPoint, expected) );

  // Touching second segment vertex
  secondSegmentPoint = Vector3d::Constant ( 0.0 );
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::TOUCHING );
  validate ( math::equals(intersectionPoint, expected) );

  // Touching first segment vertex
  firstSegmentPoint = Vector3d::Constant ( 0.0 );
  secondSegmentPoint = Vector3d::Constant ( -1.0 );
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::TOUCHING );
  validate ( math::equals(intersectionPoint, expected) );

  // Parallel segment with distance to plain
  firstSegmentPoint << 0.0, 0.0, -3.0;
  intersectionPoint << 1.0, 2.0, 3.0; // should not be modified
  expected = intersectionPoint;
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::NO_INTERSECTION );
  validate ( math::equals(intersectionPoint, expected) );

  // Parallel segment contained in plane
  firstSegmentPoint << 0.0, 0.0, 0.0;
  secondSegmentPoint << 1.0, 1.0, -2.0;
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::CONTAINED );
  validate ( math::equals(intersectionPoint, expected) );

  // Segment ending before intersection
  firstSegmentPoint << -2.0, -2.0, -2.0;
  secondSegmentPoint << -1.0, -1.0, -1.0;
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::NO_INTERSECTION );
  validate ( math::equals(intersectionPoint, expected) );

  // Segment ending before intersection (inversed segment points)
  firstSegmentPoint <<-1.0, -1.0, -1.0;
  secondSegmentPoint << -2.0, -2.0, -2.0;
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::NO_INTERSECTION );
  validate ( math::equals(intersectionPoint, expected) );
}

void GeometryComputationsTest:: testProjectVector ()
{
  TRACE();
  Eigen::Vector3d vector3D ( 1.0, 2.0, 3.0 );
  Eigen::Vector2d vector2D;
  Eigen::Vector2d vectorExpected ( 1.0, 2.0 );

  vector2D = GeometryComputations::projectVector ( vector3D, 2 );
  validate ( math::equals(vector2D, vectorExpected) );

  vector2D = GeometryComputations::projectVector ( vector3D, 1 );
  vectorExpected << 1.0, 3.0;
  validate ( math::equals(vector2D, vectorExpected) );

  vector2D = GeometryComputations::projectVector ( vector3D, 0 );
  vectorExpected << 2.0, 3.0;
  validate ( math::equals(vector2D, vectorExpected) );
}

void GeometryComputationsTest:: testContainedInTriangle ()
{
  TRACE();
  Eigen::Vector2d triangleVertex0 ( 0.0, 0.0 );
  Eigen::Vector2d triangleVertex1 ( 1.0, 0.0 );
  Eigen::Vector2d triangleVertex2 ( 0.0, 1.0 );
  Eigen::Vector2d point ( 0.25, 0.25 );

  // Contained point
  int result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::CONTAINED );

  // Not contained points
  point << -1.0, -1.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );
  point << 1.0, 1.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );
  point << 2.0, 0.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );
  point << 0.0, 2.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  // Touching points
  point << 0.0, 0.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  point << 1.0, 0.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  point << 0.0, 1.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  point << 0.5, 0.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  point << 0.0, 0.5;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  point << 0.5, 0.5;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
}

void GeometryComputationsTest:: testContainedInHyperrectangle ()
{
  TRACE();
  // 2D
  Eigen::Vector2d center2D(0, 0);
  Eigen::Vector2d sidelengths2D(1, 1);
  
  // Not contained 2D
  Eigen::Vector2d testPoint2D(2, 2);
  int result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint2D << -2.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint2D << 2.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint2D << -2.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint2D << 2.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint2D << -2.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint2D << 0.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint2D << 0.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  // Contained 2D
  testPoint2D << 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint2D << 0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint2D << -0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint2D << 0.25, -0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint2D << -0.25, -0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint2D << 0.49999999999, 0.49999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint2D << -0.49999999999, -0.49999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  // Touching 2D
  testPoint2D << 0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << -0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.5, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << -0.5, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.4999999999999999, 0.4999999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.500000000000001, 0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << -0.500000000000001, -0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << -0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.0, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.0, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << -0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.0, 0.499999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.0, -0.499999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.0, 0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.0, -0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.499999999999999, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << -0.499999999999999, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << 0.500000000000001, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint2D << -0.500000000000001, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  // 3D
  Eigen::Vector3d center3D = Eigen::Vector3d::Zero();
  Eigen::Vector3d sidelengths3D = Eigen::Vector3d::Constant(1.0);

  // Not contained 3D
  Eigen::Vector3d testPoint3D(2, 2, 2);
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint3D << -2.0, -2.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint3D << 2.0, -2.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint3D << -2.0, 2.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint3D << 2.0, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint3D << -2.0, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint3D << 0.0, 2.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint3D << 0.0, 0.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  testPoint3D << 0.0, 0.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  // Contained 3D
  testPoint3D << 0.0, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint3D << 0.25, 0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint3D << -0.25, 0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint3D << 0.25, -0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint3D << -0.25, -0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint3D << 0.25, 0.25, -0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint3D << 0.49999999999, 0.49999999999, 0.49999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  testPoint3D << -0.49999999999, -0.49999999999, -0.49999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  // Touching 3D
  testPoint3D << 0.5, 0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << -0.5, 0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.5, -0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << -0.5, -0.5, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.4999999999999999, 0.4999999999999999, 0.4999999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.500000000000001, 0.500000000000001, 0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << -0.500000000000001, -0.500000000000001, -0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.5, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << -0.5, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, 0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, -0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.5, 0.0, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << -0.5, 0.0, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, 0.499999999999999, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, -0.499999999999999, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, 0.500000000000001, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, -0.500000000000001, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.499999999999999, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << -0.499999999999999, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.500000000000001, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << -0.500000000000001, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, 0.0, 0.499999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, 0.0, -0.499999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, 0.0, 0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  testPoint3D << 0.0, 0.0, -0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );
}

}}} // namespace precice, utils, tests

