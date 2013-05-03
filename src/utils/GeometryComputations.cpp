// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "GeometryComputations.hpp"
#include "tarch/la/VectorVectorOperations.h"

namespace precice {
namespace utils {

//bool GeometryComputations:: segmentsIntersect
//(
//  const tarch::la::Vector<2,double> & a,
//  const tarch::la::Vector<2,double> & b,
//  const tarch::la::Vector<2,double> & c,
//  const tarch::la::Vector<2,double> & d,
//  bool countTouchingAsIntersection )
//{
////  precicePrint ( "segmentsIntersect: a = " << a << ", b = " << b << ", c = "
////                 << c << ", d = " << d );
//  if ( countTouchingAsIntersection ) {
//    if ( between(a, b, c) ) {
//      return true;
//    }
//    else if ( between(a, b, d) ) {
//      return true;
//    }
//    else if ( between(c, d, a) ) {
//      return true;
//    }
//    else if ( between(c, d, b) ) {
//      return true;
//    }
//  }
//
//  double abc = triangleArea(a, b, c);
//  double abd = triangleArea(a, b, d);
//  double cda = triangleArea(c, d, a);
//  double cdb = triangleArea(c, d, b);
//
//  using tarch::la::norm2;
//  double circABC = norm2(a-b) + norm2(b-c) + norm2(c-a);
//  double circABD = norm2(a-b) + norm2(b-d) + norm2(d-a);
//  double circCDA = norm2(c-d) + norm2(d-a) + norm2(a-c);
//  double circCDB = norm2(c-d) + norm2(d-b) + norm2(b-c);
//
////  precicePrint ( "triA(abc) = " << abc << ", triA(abd) = " << abd
////                 << ", triA(cda) = " << cda << ", triA(cdb) = " << cdb );
////  precicePrint ( "circ(abc) = " << circABC
////                 << ", circ(abd) = " << circABD
////                 << ", circ(cda) = " << circCDA
////                 << ", circ(cdb) = " << circCDB );
//
//  abc /= circABC;
//  abd /= circABD;
//  cda /= circCDA;
//  cdb /= circCDB;
//
////  precicePrint ( "after scaling: triA(abc) = " << abc << ", triA(abd) = " << abd
////                 << ", triA(cda) = " << cda << ", triA(cdb) = " << cdb );
//
//  // Check, if one point lies on line defined by segment and the other not (-> xor).
//  // If true, either one point is between, which means the segments are
//  // touching only, or the segments are neither touching nor intersecting
//  // (-> false). This case of touching segments is detected in the beginning
//  // of this function, if countTouchingAsIntersection is true. Otherwise,
//  // it should not be counted (-> false).
//  if ( xOR(std::abs(abc) <= tarch::la::NUMERICAL_ZERO_DIFFERENCE,
//       std::abs(abd) <= tarch::la::NUMERICAL_ZERO_DIFFERENCE) )
//  {
////    precicePrint ( "return at 1" );
//    return false;
//  }
//  if ( xOR(std::abs(cda) <= tarch::la::NUMERICAL_ZERO_DIFFERENCE,
//       std::abs(cdb) <= tarch::la::NUMERICAL_ZERO_DIFFERENCE) )
//  {
////    precicePrint ( "return at 2" );
//    return false;
//  }
//
//  // Check, whether the segments are intersecting in the real sense.
//  bool isFirstSegmentBetween = xOR (abc > - tarch::la::NUMERICAL_ZERO_DIFFERENCE,
//                                    abd > - tarch::la::NUMERICAL_ZERO_DIFFERENCE );
//  bool isSecondSegmentBetween = xOR (cda > - tarch::la::NUMERICAL_ZERO_DIFFERENCE,
//                                     cdb > - tarch::la::NUMERICAL_ZERO_DIFFERENCE );
////  precicePrint ( "return at 3 = " << isFirstSegmentBetween << " && "
////                 << isSecondSegmentBetween );
//  return isFirstSegmentBetween && isSecondSegmentBetween;
//}

bool GeometryComputations:: lineIntersection
(
   const tarch::la::Vector<2,double> & a,
   const tarch::la::Vector<2,double> & b,
   const tarch::la::Vector<2,double> & c,
   const tarch::la::Vector<2,double> & d,
   tarch::la::Vector<2,double>       & intersectionPoint )
{
   // Compute denominator for solving 2x2 equation system
   double D = a(0)*(d(1)-c(1)) +
              b(0)*(c(1)-d(1)) +
              d(0)*(b(1)-a(1)) -
              c(0)*(a(1)-b(1));

   // If D==0, the two lines are parallel
   if ( tarch::la::equals(D, 0.0) ) {
      return false;
   }

   // Compute line segment parameter s, which is in [0, 1], if the
   // intersection point is within the segment.
   double s = (a(0)*(d(1)-c(1)) + c(0)*(a(1)-d(1)) + d(0)*(c(1)-a(1))) / D;

   intersectionPoint = a + s * (b - a);
   return true;
}

int GeometryComputations:: segmentPlaneIntersection
(
   const tarch::la::Vector<3,double> & pointOnPlane,
   const tarch::la::Vector<3,double> & planeNormal,
   const tarch::la::Vector<3,double> & firstPointSegment,
   const tarch::la::Vector<3,double> & secondPointSegment,
   tarch::la::Vector<3,double>       & intersectionPoint )
{
   // Methodology of "Computation Geometry in C", Joseph O'Rourke, Chapter 7.3.1

   // The plane is represented by: p(x,y,z) dot planeNormal = d  (1)
   // We need to compute d by using pointOnPlane as p(x,y,z):
   double d = pointOnPlane * planeNormal;

   // The segment is represented by:
   // firstPointSegment + (secondPointSegment - firstPointSegement) * t (2)
   // We insert (2) into (1) and reorder for t. The resulting equations reads
   // t = (d - firstPointSegment dot planeNormal) /
   //     ((secondPointSegment - firstPointSegment) dot planeNormal). (3)
   //
   // If the denominator of (3) equals 0, the segment is parallel to the plane.
   // If, in addition, the denominator of (3) equals 0, the segment lies in the
   // plane. Otherwise we can compute a point of intersection.
   tarch::la::Vector<3,double> segmentVec ( secondPointSegment - firstPointSegment );
   double nominator = d - firstPointSegment * planeNormal;
   double denominator = segmentVec * planeNormal;
   if ( tarch::la::equals(denominator, 0.0) ) {
      if ( tarch::la::equals(nominator, 0.0) ) {
         return CONTAINED;
      }
      else {
         return NO_INTERSECTION;
      }
   }
   double t = nominator / denominator;

   // If t is larger than 1 or smaller than zero, the intersection is not within
   // the line segment
   if ( tarch::la::greater(t, 1.0) || tarch::la::greater(0.0, t) ) {
      return NO_INTERSECTION;
   }

   intersectionPoint = firstPointSegment + segmentVec * t;

   // If t equals 1 or 0, the segment is just touching the plane, otherwise, a
   // real intersection is happening.
   if ( tarch::la::equals(t, 0.0) || tarch::la::equals(t, 1.0) ) {
      return TOUCHING;
   }

   return INTERSECTION;
}

//double GeometryComputations:: triangleArea
//(
//   const tarch::la::Vector<2,double> & a,
//   const tarch::la::Vector<2,double> & b,
//   const tarch::la::Vector<2,double> & c )
//{
//   tarch::la::Vector<2,double> A = b - a;
//   tarch::la::Vector<2,double> B = c - a;
//   return 0.5 * (A(0)*B(1) - A(1)*B(0));
//}


//double GeometryComputations:: triangleArea
//(
//   const tarch::la::Vector<3, double> & a,
//   const tarch::la::Vector<3, double> & b,
//   const tarch::la::Vector<3, double> & c )
//{
//   tarch::la::Vector<3,double> A = b - a;
//   tarch::la::Vector<3,double> B = c - a;
//   tarch::la::Vector<3,double> result;
//
//   return 0.5 * tarch::la::norm2 ( tarch::la::cross(A,B,result) );
//}

double GeometryComputations:: tetraVolume
(
  const tarch::la::Vector<3,double> & a,
  const tarch::la::Vector<3,double> & b,
  const tarch::la::Vector<3,double> & c,
  const tarch::la::Vector<3,double> & d )
{
  // TODO compute volume of tetraheder
  return -1.0;
}

tarch::la::Vector<2,double> GeometryComputations:: projectVector
(
   const tarch::la::Vector<3,double> & vector,
   int indexDimensionToRemove )
{
   assertion ( indexDimensionToRemove >= 0 );
   assertion ( indexDimensionToRemove < 3 );
   tarch::la::Vector<2, double> projectedVector;
   int reducedDim = 0;
   for (int dim=0; dim < 3; dim++) {
      if ( indexDimensionToRemove == dim ) {
         continue;
      }
      projectedVector(reducedDim) = vector(dim);
      reducedDim++;
   }
   return projectedVector;
}

int GeometryComputations:: containedInTriangle
(
  const tarch::la::Vector<2,double> & triangleVertex0,
  const tarch::la::Vector<2,double> & triangleVertex1,
  const tarch::la::Vector<2,double> & triangleVertex2,
  const tarch::la::Vector<2,double> & testPoint )
{
  using namespace tarch::la;
  double area0 = triangleArea ( triangleVertex0, triangleVertex1, testPoint );
  double area1 = triangleArea ( triangleVertex1, triangleVertex2, testPoint );
  double area2 = triangleArea ( triangleVertex2, triangleVertex0, testPoint );

  double length01 = norm2 ( triangleVertex0 - triangleVertex1 );
  double length12 = norm2 ( triangleVertex1 - triangleVertex2 );
  double length20 = norm2 ( triangleVertex2 - triangleVertex0 );
  double length1t = norm2 ( triangleVertex1 - testPoint );
  double length2t = norm2 ( triangleVertex2 - testPoint );
  double length0t = norm2 ( triangleVertex0 - testPoint );
  double scale0 = length01 + length1t + length0t;
  double scale1 = length12 + length2t + length1t;
  double scale2 = length20 + length0t + length2t;
  area0 /= scale0;
  area1 /= scale1;
  area2 /= scale2;

  int sign0 = sign(area0);
  int sign1 = sign(area1);
  int sign2 = sign(area2);
  int absSumSigns = std::abs( sign0 + sign1 + sign2 );
  if ( absSumSigns == 3 ) {
    // In triangle
    return CONTAINED;
  }
  else if ( absSumSigns == 2 ) {
    // On triangle edge
    return TOUCHING;
  }
  else if ( std::abs(sign0) + std::abs(sign1) + std::abs(sign2) == 1 ) {
    // On triangle vertex
    return TOUCHING;
  }
  // Outside of triangle
  return NOT_CONTAINED;
}

}} // namespace precice, utils
