// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "GeometryComputationsTest.hpp"
#include "../GeometryComputations.hpp"
#include "../Dimensions.hpp"
#include "../Parallel.hpp"
#include "../Helpers.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::GeometryComputationsTest)

namespace precice {
namespace utils {
namespace tests {

tarch::logging::Log GeometryComputationsTest::
    _log ("precice::utils::GeometryComputationsTest");

GeometryComputationsTest:: GeometryComputationsTest ()
:
  TestCase ("utils::GeometryComputationsTest")
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
  preciceTrace ( "testCollinear()" );
  // 2D test setup
  tarch::la::Vector<2,double> a2D (0.0);
  tarch::la::Vector<2,double> b2D (1.0);
  tarch::la::Vector<2,double> collinearPoint2D (0.5);
  tarch::la::Vector<2,double> notCollinearPoint2D (0.5, 0.6);

  // 3D test setup
  tarch::la::Vector<3,double> a3D (0.0);
  tarch::la::Vector<3,double> b3D (1.0);
  tarch::la::Vector<3,double> collinearPoint3D (0.5);
  tarch::la::Vector<3,double> notCollinearPoint3D (0.5, 0.6, 0.5);

  // 2D test validations
  validate ( GeometryComputations::collinear(a2D, b2D, collinearPoint2D) );
  validate ( ! GeometryComputations::collinear(a2D, b2D, notCollinearPoint2D) );

  // 3D test validations
  validate ( GeometryComputations::collinear(a3D, b3D, collinearPoint3D) );
  validate ( ! GeometryComputations::collinear(a3D, b3D, notCollinearPoint3D) );
}

void GeometryComputationsTest:: testBetween ()
{
  preciceTrace ( "testBetween()" );
  for ( int dim=2; dim <= 3; dim++ ){
    utils::DynVector a(dim);
    utils::DynVector b(dim);
    utils::DynVector betweenPoint(dim);
    utils::DynVector betweenLimitPoint(dim);
    utils::DynVector collinearOutsidePoint(dim);
    utils::DynVector outsidePoint(dim);
    if ( dim == 2 ){
      assignList(a) = 0.0, 0.0;
      assignList(b) = 1.0, 1.0;
      assignList(betweenPoint) = 0.5, 0.5;
      assignList(betweenLimitPoint) = 1.0, 1.0;
      assignList(collinearOutsidePoint) = 2.0, 2.0;
      assignList(outsidePoint) = 0.5, 0.4;
    }
    else {
      assignList(a) = 0.0, 0.0, 0.0;
      assignList(b) = 1.0, 1.0, 1.0;
      assignList(betweenPoint) = 0.5, 0.5, 0.5;
      assignList(betweenLimitPoint) = 1.0, 1.0, 1.0;
      assignList(collinearOutsidePoint) = 2.0, 2.0, 2.0;
      assignList(outsidePoint) = 0.5, 0.4, 0.5;
    }
    validate(utils::GeometryComputations::between(a, b, betweenPoint));
    validate(utils::GeometryComputations::between(a, b, betweenLimitPoint));
    validate(! utils::GeometryComputations::between(a, b, collinearOutsidePoint));
    validate(! utils::GeometryComputations::between(a, b, outsidePoint));
  }
}

void GeometryComputationsTest:: testTriangleArea ()
{
  preciceTrace ( "testTriangleArea()" );
  { // 2D
    utils::Vector2D a, b, c;
    double area;
    assignList(a) = 0.0, 0.0;
    assignList(b) = 1.0, 0.0;
    assignList(c) = 0.0, 1.0;
    area = utils::GeometryComputations::triangleArea(a, b, c);
    validateNumericalEquals (area, 0.5);

    assignList(b) = 0.0, 1.0;
    assignList(c) = 1.0, 0.0;
    area = utils::GeometryComputations::triangleArea(a, b, c);
    validateNumericalEquals (area, -0.5);
  }
  { // 3D
    utils::Vector3D a, b, c;
    double area;
    assignList(a) = 0.0, 0.0, 0.0;
    assignList(b) = 1.0, 0.0, 1.0;
    assignList(c) = 1.0, 1.0, 1.0;
    area = utils::GeometryComputations::triangleArea(a, b, c);
    validateNumericalEquals (area, 0.5 * sqrt(2.0));
  }
}

void GeometryComputationsTest:: testSegmentsIntersect ()
{
  preciceTrace ( "testSegmentsIntersect" );
  typedef GeometryComputations GeoComp;
  utils::Vector2D a(0.0), b(0.0), c(0.0), d(0.5);
  assignList(b) = 1.0, 0.0;
  assignList(c) = 0.5, 0.0;
  validate ( GeoComp::segmentsIntersect(a,b,c,d,true));
  validate ( ! GeoComp::segmentsIntersect(a,b,c,d,false));

  c(1) = -0.2;
  validate ( GeoComp::segmentsIntersect(a,b,c,d,false));

  // Test case motivated from bug
  assignList(a) = 0.23104429, 1.87753905;
  assignList(b) = 0.22608634, 1.88114120;
  assignList(c) = 0.23058985, 1.87733882;
  assignList(d) = 0.23058985, 1.87846349;
  validate ( GeoComp::segmentsIntersect(a,b,c,d,false));

  // Another bug motivated test, T-intersection slightly above numerical accuracy
  assignList(a) = 1.0, 1.0;
  assignList(b) = 0.0, 1.0;
  assignList(c) = 1.0 - (10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE), 1.1;
  assignList(d) = 1.0 - (10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE), 0.9;
  validate ( GeoComp::segmentsIntersect(a,b,c,d,false) );
  validate ( GeoComp::segmentsIntersect(c,d,a,b,false) );
}

void GeometryComputationsTest:: testSegmentPlaneIntersection ()
{
  preciceTrace ( "testSegmentPlaneIntersection()" );
  using utils::Vector3D;
  Vector3D planeNormal ( 1.0 );
  Vector3D pointOnPlane ( 0.0 );
  Vector3D firstSegmentPoint ( 1.0 );
  Vector3D secondSegmentPoint ( -1.0 );
  Vector3D intersectionPoint ( 1.0 );
  Vector3D expected ( 0.0 );

  // True intersection
  int result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::INTERSECTION );
  validate ( tarch::la::equals(intersectionPoint, expected) );

  // Touching second segment vertex
  secondSegmentPoint = Vector3D ( 0.0 );
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::TOUCHING );
  validate ( tarch::la::equals(intersectionPoint, expected) );

  // Touching first segment vertex
  firstSegmentPoint = Vector3D ( 0.0 );
  secondSegmentPoint = Vector3D ( -1.0 );
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::TOUCHING );
  validate ( tarch::la::equals(intersectionPoint, expected) );

  // Parallel segment with distance to plain
  firstSegmentPoint = Vector3D ( 0.0, 0.0, -3.0 );
  intersectionPoint = Vector3D ( 1.0, 2.0, 3.0 ); // should not be modified
  expected = intersectionPoint;
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::NO_INTERSECTION );
  validate ( tarch::la::equals(intersectionPoint, expected) );

  // Parallel segment contained in plane
  firstSegmentPoint = Vector3D ( 0.0, 0.0, 0.0 );
  secondSegmentPoint = Vector3D ( 1.0, 1.0, -2.0 );
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::CONTAINED );
  validate ( tarch::la::equals(intersectionPoint, expected) );

  // Segment ending before intersection
  firstSegmentPoint = Vector3D ( -2.0, -2.0, -2.0 );
  secondSegmentPoint = Vector3D ( -1.0, -1.0, -1.0 );
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::NO_INTERSECTION );
  validate ( tarch::la::equals(intersectionPoint, expected) );

  // Segment ending before intersection (inversed segment points)
  firstSegmentPoint = Vector3D ( -1.0, -1.0, -1.0 );
  secondSegmentPoint = Vector3D ( -2.0, -2.0, -2.0 );
  result = GeometryComputations::segmentPlaneIntersection (
      pointOnPlane, planeNormal, firstSegmentPoint,
      secondSegmentPoint, intersectionPoint );
  validateEquals ( result, GeometryComputations::NO_INTERSECTION );
  validate ( tarch::la::equals(intersectionPoint, expected) );
}

void GeometryComputationsTest:: testProjectVector ()
{
  preciceTrace ( "testProjectVector()" );
  tarch::la::Vector<3,double> vector3D ( 1.0, 2.0, 3.0 );
  tarch::la::Vector<2,double> vector2D;
  tarch::la::Vector<2,double> vectorExpected ( 1.0, 2.0 );

  vector2D = GeometryComputations::projectVector ( vector3D, 2 );
  validate ( tarch::la::equals(vector2D, vectorExpected) );

  vector2D = GeometryComputations::projectVector ( vector3D, 1 );
  assignList(vectorExpected) = 1.0, 3.0;
  validate ( tarch::la::equals(vector2D, vectorExpected) );

  vector2D = GeometryComputations::projectVector ( vector3D, 0 );
  assignList(vectorExpected) = 2.0, 3.0;
  validate ( tarch::la::equals(vector2D, vectorExpected) );
}

void GeometryComputationsTest:: testContainedInTriangle ()
{
  preciceTrace ( "testContainedInTriangle()" );
  tarch::la::Vector<2,double> triangleVertex0 ( 0.0, 0.0 );
  tarch::la::Vector<2,double> triangleVertex1 ( 1.0, 0.0 );
  tarch::la::Vector<2,double> triangleVertex2 ( 0.0, 1.0 );
  tarch::la::Vector<2,double> point ( 0.25, 0.25 );

  // Contained point
  int result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::CONTAINED );

  // Not contained points
  assignList(point) = -1.0, -1.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );
  assignList(point) = 1.0, 1.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );
  assignList(point) = 2.0, 0.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );
  assignList(point) = 0.0, 2.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  // Touching points
  assignList(point) = 0.0, 0.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  assignList(point) = 1.0, 0.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  assignList(point) = 0.0, 1.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  assignList(point) = 0.5, 0.0;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  assignList(point) = 0.0, 0.5;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
  assignList(point) = 0.5, 0.5;
  result = GeometryComputations::containedInTriangle (
      triangleVertex0, triangleVertex1, triangleVertex2, point );
  validateEquals ( result, GeometryComputations::TOUCHING );
}

void GeometryComputationsTest:: testContainedInHyperrectangle ()
{
  preciceTrace ( "testContainedInHyperrectangle()" );
  // 2D
  tarch::la::Vector<2, double> center2D ( 0.0, 0.0 );
  tarch::la::Vector<2, double> sidelengths2D ( 1.0, 1.0 );

  // Not contained 2D
  tarch::la::Vector<2, double> testPoint2D ( 2.0, 2.0 );
  int result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint2D) = -2.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint2D) = 2.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint2D) = -2.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint2D) = 2.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint2D) = -2.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint2D) = 0.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint2D) = 0.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  // Contained 2D
  assignList(testPoint2D) = 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint2D) = 0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint2D) = -0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint2D) = 0.25, -0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint2D) = -0.25, -0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint2D) = 0.49999999999, 0.49999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint2D) = -0.49999999999, -0.49999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  // Touching 2D
  assignList(testPoint2D) = 0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = -0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.5, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = -0.5, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.4999999999999999, 0.4999999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.500000000000001, 0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = -0.500000000000001, -0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = -0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.0, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.0, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = -0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.0, 0.499999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.0, -0.499999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.0, 0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.0, -0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.499999999999999, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = -0.499999999999999, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = 0.500000000000001, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint2D) = -0.500000000000001, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths2D, center2D, testPoint2D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  // 3D
  tarch::la::Vector<3, double> center3D ( 0.0 );
  tarch::la::Vector<3, double> sidelengths3D ( 1.0 );

  // Not contained 3D
  tarch::la::Vector<3, double> testPoint3D ( 2.0, 2.0, 2.0 );
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint3D) = -2.0, -2.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint3D) = 2.0, -2.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint3D) = -2.0, 2.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint3D) = 2.0, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint3D) = -2.0, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint3D) = 0.0, 2.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint3D) = 0.0, 0.0, 2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  assignList(testPoint3D) = 0.0, 0.0, -2.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::NOT_CONTAINED );

  // Contained 3D
  assignList(testPoint3D) = 0.0, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint3D) = 0.25, 0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint3D) = -0.25, 0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint3D) = 0.25, -0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint3D) = -0.25, -0.25, 0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint3D) = 0.25, 0.25, -0.25;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint3D) = 0.49999999999, 0.49999999999, 0.49999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  assignList(testPoint3D) = -0.49999999999, -0.49999999999, -0.49999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::CONTAINED );

  // Touching 3D
  assignList(testPoint3D) = 0.5, 0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = -0.5, 0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.5, -0.5, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = -0.5, -0.5, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.4999999999999999, 0.4999999999999999, 0.4999999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.500000000000001, 0.500000000000001, 0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = -0.500000000000001, -0.500000000000001, -0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.5, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = -0.5, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, 0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, -0.5, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.5, 0.0, -0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = -0.5, 0.0, 0.5;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, 0.499999999999999, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, -0.499999999999999, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, 0.500000000000001, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, -0.500000000000001, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.499999999999999, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = -0.499999999999999, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.500000000000001, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = -0.500000000000001, 0.0, 0.0;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, 0.0, 0.499999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, 0.0, -0.499999999999999;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, 0.0, 0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );

  assignList(testPoint3D) = 0.0, 0.0, -0.500000000000001;
  result = GeometryComputations::containedInHyperrectangle (
    sidelengths3D, center3D, testPoint3D );
  validateEquals ( result, GeometryComputations::TOUCHING );
}

}}} // namespace precice, utils, tests

