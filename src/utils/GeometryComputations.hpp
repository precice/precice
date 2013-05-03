// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_GEOMETRYCOMPUTATIONS_HPP_
#define PRECICE_UTILS_GEOMETRYCOMPUTATIONS_HPP_

#include "Dimensions.hpp"
#include "Globals.hpp"
#include "boost/static_assert.hpp"
#include "tarch/la/VectorOperations.h"
#include "tarch/utils/EnableIf.h"
#include "utils/Helpers.hpp"
#include <iostream>

namespace precice {
namespace utils {


/**
 * @brief Provides computational geometry operations.
 */
class GeometryComputations
{
public:

   enum ResultConstants {
      NO_INTERSECTION,
      INTERSECTION,
      TOUCHING,
      NOT_CONTAINED,
      CONTAINED,
   };

   /**
    * @brief Determines, if to line segments intersect properly or inproperly
    *
    * Works only for Dim2.
    * The method first checks, if one of the segment end points lies on the
    * other segment. If yes, this is inproper intersection, and considered
    * to be valid.
    * Then, it checks, if it is true for both segments, that one point of
    * the other segment lies right of the segment and the other point left
    * of it. Then, the segments do intersect each other properly.
    *
    * @param a [IN] First point of line segment 1
    * @param b [IN] Second point of line segement 1
    * @param c [IN] First point of line segment 2
    * @param d [IN] Second point of line segment 2
    *
    * @return True, if interseting. False, otherwise
    */
   template<typename VECTORA_T, typename VECTORB_T, typename VECTORC_T, typename VECTORD_T>
   static bool segmentsIntersect (
      const VECTORA_T& a,
      const VECTORB_T& b,
      const VECTORC_T& c,
      const VECTORD_T& d,
      bool countTouchingAsIntersection );

   /**
    * @brief Determines the intersection point of two lines, if one exists.
    *
    * Works only for Dim2.
    * @param a [IN] First point on first line.
    * @param b [IN] Second point on first line.
    * @param c [IN] First point on second line.
    * @param d [IN] Second point on second line.
    * @param intersectionPoint [OUT] Point of intersection.
    * @return True, if an intersection point exists and could be determined.
    */
   static bool lineIntersection (
      const tarch::la::Vector<2,double> & a,
      const tarch::la::Vector<2,double> & b,
      const tarch::la::Vector<2,double> & c,
      const tarch::la::Vector<2,double> & d,
      tarch::la::Vector<2,double>       & intersectionPoint );

   /**
    * @brief Determines the intersection point of a segment with a plane in 3D.
    *
    * @param intersectionPoint [OUT] Contains coordinates of intersection point,
    *        if return value equals TOUCHING or INTERSECTION.
    *
    * @return Returns type of intersection. One of
    *         - NO_INTERSECTION
    *         - INTERSECTION
    *         - TOUCHING
    *         - CONTAINED
    */
   static int segmentPlaneIntersection (
      const tarch::la::Vector<3,double> & pointOnPlane,
      const tarch::la::Vector<3,double> & planeNormal,
      const tarch::la::Vector<3,double> & firstPointSegment,
      const tarch::la::Vector<3,double> & secondPointSegment,
      tarch::la::Vector<3,double>       & intersectionPoint );

   /**
    * @brief Determines, if a point lies on the line segment defined by a, b
    *
    * Works for 2D and 3D.
    *
    * @param a       [IN] First point of line segment
    * @param b       [IN] Second point of line segment
    * @param toCheck [IN] Point to be checked for "betweenness"
    *
    * @return True, if toCheck lies between a and b. False, otherwise
    */
   template<typename VECTORA_T, typename VECTORB_T, typename VECTORC_T>
   static bool between (
      const VECTORA_T& a,
      const VECTORB_T& b,
      const VECTORC_T& toCheck );

   /**
    * @brief Determines, if three points are collinear (on one line)
    *
    * Works for Dim2 and Dim3.
    *
    * @param a [IN] First point to check
    * @param b [IN] Second point to check
    * @param c [IN] Third point to check
    *
    * @return True, if collinear. False, otherwise
    */
   template<typename VECTORA_T, typename VECTORB_T, typename VECTORC_T>
   static bool collinear (
      const VECTORA_T& a,
      const VECTORB_T& b,
      const VECTORC_T& c );

   /**
    * @brief Determines, if two lines are parallel to each other.
    *
    * Works for Dim2 and Dim3.
    *
    * @param a [IN] First point on first line.
    * @param b [IN] Second point on first line.
    * @param c [IN] First point on second line.
    * @param d [IN] Second point on second line.
    * @return True, if the two lines are parallel to each other.
    */
   template< int dim >
   static bool parallel (
      const tarch::la::Vector<dim,double> & a,
      const tarch::la::Vector<dim,double> & b,
      const tarch::la::Vector<dim,double> & c,
      const tarch::la::Vector<dim,double> & d );

   /**
    * @brief Computes the signed area of a triangle in 2D.
    *
    * The area is negative, if the points a, b, c are given in
    * clockwise ordering, otherwise positive.
    */
   template<typename VECTOR>
   static typename tarch::utils::EnableIf< tarch::la::IsVector<VECTOR>::value,
     double
   >::Type triangleArea (
      const VECTOR& a,
      const VECTOR& b,
      const VECTOR& c );

   /**
    * @brief Computes the (unsigned) area of a triangle in 3D.
    */
//   static double triangleArea (
//     const tarch::la::Vector<3,double> & a,
//     const tarch::la::Vector<3,double> & b,
//     const tarch::la::Vector<3,double> & c );

   static double tetraVolume (
     const tarch::la::Vector<3,double> & a,
     const tarch::la::Vector<3,double> & b,
     const tarch::la::Vector<3,double> & c,
     const tarch::la::Vector<3,double> & d );

   /**
    * @brief Projects a 3D vector to a 2D one by removing one dimension.
    */
   static tarch::la::Vector<2, double> projectVector (
      const tarch::la::Vector<3, double> & vector,
      int indexDimensionToRemove );

   /**
    * @brief Tests if a vertex is contained in a triangle.
    *
    * @return One of:
    *         - CONTAINED
    *         - NOT_CONTAINED
    *         - TOUCHING
    */
   static int containedInTriangle (
      const tarch::la::Vector<2,double> & triangleVertex0,
      const tarch::la::Vector<2,double> & triangleVertex1,
      const tarch::la::Vector<2,double> & triangleVertex2,
      const tarch::la::Vector<2,double> & testPoint );

   /**
    * @brief Tests, if a vertex is contained in a hyperrectangle.
    *
    * @return One of:
    *         - CONTAINED
    *         - NOT_CONTAINED
    *         - TOUCHING
    */
   template<typename VECTORA_T, typename VECTORB_T, typename VECTORC_T>
   static int containedInHyperrectangle (
      const VECTORA_T& sidelengths,
      const VECTORB_T& center,
      const VECTORC_T& testPoint );
};


// --------------------------------------------------------- HEADER DEFINITIONS

template<typename VECTORA_T, typename VECTORB_T, typename VECTORC_T, typename VECTORD_T>
bool GeometryComputations:: segmentsIntersect
(
  const VECTORA_T& a,
  const VECTORB_T& b,
  const VECTORC_T& c,
  const VECTORD_T& d,
  bool countTouchingAsIntersection )
{
  assertion1 ( a.size() == 2, a.size() );
  assertion1 ( b.size() == 2, b.size() );
  assertion1 ( c.size() == 2, c.size() );
  assertion1 ( d.size() == 2, d.size() );

  if ( countTouchingAsIntersection ) {
    if ( between(a, b, c) ) {
      return true;
    }
    else if ( between(a, b, d) ) {
      return true;
    }
    else if ( between(c, d, a) ) {
      return true;
    }
    else if ( between(c, d, b) ) {
      return true;
    }
  }

  double abc = triangleArea(a, b, c);
  double abd = triangleArea(a, b, d);
  double cda = triangleArea(c, d, a);
  double cdb = triangleArea(c, d, b);

  using tarch::la::norm2;
  double circABC = norm2(a-b) + norm2(b-c) + norm2(c-a);
  double circABD = norm2(a-b) + norm2(b-d) + norm2(d-a);
  double circCDA = norm2(c-d) + norm2(d-a) + norm2(a-c);
  double circCDB = norm2(c-d) + norm2(d-b) + norm2(b-c);

  abc /= circABC;
  abd /= circABD;
  cda /= circCDA;
  cdb /= circCDB;

  // Check, if one point lies on line defined by segment and the other not (-> xor).
  // If true, either one point is between, which means the segments are
  // touching only, or the segments are neither touching nor intersecting
  // (-> false). This case of touching segments is detected in the beginning
  // of this function, if countTouchingAsIntersection is true. Otherwise,
  // it should not be counted (-> false).
  if ( xOR(std::abs(abc) <= tarch::la::NUMERICAL_ZERO_DIFFERENCE,
           std::abs(abd) <= tarch::la::NUMERICAL_ZERO_DIFFERENCE) )
  {
    return false;
  }
  if ( xOR(std::abs(cda) <= tarch::la::NUMERICAL_ZERO_DIFFERENCE,
           std::abs(cdb) <= tarch::la::NUMERICAL_ZERO_DIFFERENCE) )
  {
    return false;
  }

  // Check, whether the segments are intersecting in the real sense.
  bool isFirstSegmentBetween = xOR (abc > - tarch::la::NUMERICAL_ZERO_DIFFERENCE,
                                    abd > - tarch::la::NUMERICAL_ZERO_DIFFERENCE );
  bool isSecondSegmentBetween = xOR (cda > - tarch::la::NUMERICAL_ZERO_DIFFERENCE,
                                     cdb > - tarch::la::NUMERICAL_ZERO_DIFFERENCE );
  return isFirstSegmentBetween && isSecondSegmentBetween;
}

template<typename VECTORA_T, typename VECTORB_T, typename VECTORC_T>
bool GeometryComputations:: between
(
  const VECTORA_T& a,
  const VECTORB_T& b,
  const VECTORC_T& toCheck )
{
  if (not collinear(a, b, toCheck)) {
    return false;
  }

  if (a(0) != b(0)) {
    return ( tarch::la::greaterEquals(toCheck(0), a(0)) &&
             tarch::la::greaterEquals(b(0), toCheck(0)) ) ||
             ( tarch::la::greaterEquals(a(0), toCheck(0)) &&
             tarch::la::greaterEquals(toCheck(0), b(0)) );
  }
  else {
    return ( tarch::la::greaterEquals(toCheck(1), a(1)) &&
             tarch::la::greaterEquals(b(1), toCheck(1)) ) ||
             ( tarch::la::greaterEquals(a(1), toCheck(1)) &&
             tarch::la::greaterEquals(toCheck(1), b(1)) );
  }
}

//template< int dim >
//bool GeometryComputations:: collinear
//(
//   const tarch::la::Vector<dim,double> & a,
//   const tarch::la::Vector<dim,double> & b,
//   const tarch::la::Vector<dim,double> & c )

template<typename VECTORA_T, typename VECTORB_T, typename VECTORC_T>
bool GeometryComputations:: collinear (
  const VECTORA_T& a,
  const VECTORB_T& b,
  const VECTORC_T& c )
{
  assertion2 ( a.size() == b.size(), a.size(), b.size() );
  assertion2 ( a.size() == c.size(), a.size(), c.size() );
  using namespace tarch::la;
  double triangleOutline = norm2(b-a) + norm2(c-b) + norm2(a-c);
  if ( equals(triangleArea(a, b, c) / triangleOutline, 0.0) ) {
    return true;
  }
  return false;
}

template< int dim >
bool GeometryComputations:: parallel
(
   const tarch::la::Vector<dim,double> & a,
   const tarch::la::Vector<dim,double> & b,
   const tarch::la::Vector<dim,double> & c,
   const tarch::la::Vector<dim,double> & d )
{
   if ( tarch::la::equals(triangleArea(a, b, c), 0.0)
        && tarch::la::equals(triangleArea(a, b, d), 0.0) )
   {
      return true;
   }
   return false;
}

template<typename VECTOR>
typename tarch::utils::EnableIf< tarch::la::IsVector<VECTOR>::value,
  double
>::Type GeometryComputations:: triangleArea
(
   const VECTOR& a,
   const VECTOR& b,
   const VECTOR& c )
{
  assertion2 ( a.size() == b.size(), a.size(), b.size() );
  assertion2 ( b.size() == c.size(), b.size(), c.size() );
  if ( a.size() == 2 ){
    utils::Vector2D A = b;
    A -= a;
    utils::Vector2D B = c;
    B -= a;
    return 0.5 * (A(0)*B(1) - A(1)*B(0));
  }
  else {
    assertion1 ( a.size() == 3, a.size() );
    utils::Vector3D A = b; A -= a;
    utils::Vector3D B = c; B -= a;
    utils::Vector3D result;
    return 0.5 * tarch::la::norm2( tarch::la::cross(A,B,result) );
  }
}

//template< int dim >
//int GeometryComputations:: containedInHyperrectangle
//(
//  const tarch::la::Vector<dim,double> & sidelengths,
//  const tarch::la::Vector<dim,double> & center,
//  const tarch::la::Vector<dim,double> & testPoint )

template<typename VECTORA_T, typename VECTORB_T, typename VECTORC_T>
int GeometryComputations:: containedInHyperrectangle
(
  const VECTORA_T& sidelengths,
  const VECTORB_T& center,
  const VECTORC_T& testPoint )
{
  int dim = sidelengths.size();
  assertion2 ( dim == center.size(), dim, center.size() );
  assertion2 ( dim == testPoint.size(), dim, testPoint.size() );
  utils::DynVector toCenter(testPoint);
  toCenter -= center;
  tarch::la::abs ( toCenter, toCenter );

  double diff = 0.0;
  bool touching = false;
  for ( int i=0; i < dim; i++ ) {
    diff = 0.5 * sidelengths(i) - toCenter(i);
    if ( tarch::la::greater(0.0, diff) ) {
      return NOT_CONTAINED;
    }
    if ( tarch::la::equals(diff, 0.0) ) {
      touching = true;
    }
  }
  if ( touching ) {
    return TOUCHING;
  }
  return CONTAINED;
}

}} // namespace precice, utils

#endif /* PRECICE_UTILS_GEOMETRYCOMPUTATIONS_HPP_ */
