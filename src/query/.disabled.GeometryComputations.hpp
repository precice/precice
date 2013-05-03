// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QUERY_GEOMETRYCOMPUTATIONS_HPP_
#define PRECICE_QUERY_GEOMETRYCOMPUTATIONS_HPP_

namespace precice {
namespace query {

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
template< typename Vector >
typename boost::enable_if_c< la::IsVector<Vector>::value,
bool>::type segmentsIntersect (
  const Vector & a,
  const Vector & b,
  const Vector & c,
  const Vector & d,
  bool           countTouchingAsIntersection )
{
  typedef la::VectorTraits<Vector> Traits;
  assertion ( Traits::size(a) == 2 );
  assertion ( Traits::size(b) == 2 );
  assertion ( Traits::size(c) == 2 );
  assertion ( Traits::size(d) == 2 );

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

  double circABC = la::norm(a-b) + la::norm(b-c) + la::norm(c-a);
  double circABD = la::norm(a-b) + la::norm(b-d) + la::norm(d-a);
  double circCDA = la::norm(c-d) + la::norm(d-a) + la::norm(a-c);
  double circCDB = la::norm(c-d) + la::norm(d-b) + la::norm(b-c);

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
  if ( xOR(la::equals(abc, 0.0), la::equals(abd, 0.0)) ) {
    return false;
  }
  if ( xOR(la::equals(cda, 0.0) , la::equals(cdb, 0.0)) ) {
    return false;
  }

  // Check, whether the segments are intersecting in the real sense.
  bool isFirstSegmentBetween = xOR (abc > - tarch::la::NUMERICAL_ZERO_DIFFERENCE,
                                    abd > - tarch::la::NUMERICAL_ZERO_DIFFERENCE );
  bool isSecondSegmentBetween = xOR (cda > - tarch::la::NUMERICAL_ZERO_DIFFERENCE,
                                     cdb > - tarch::la::NUMERICAL_ZERO_DIFFERENCE );
  return isFirstSegmentBetween && isSecondSegmentBetween;
}

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
bool lineIntersection (
   const tarch::la::Vector<2, double> & a,
   const tarch::la::Vector<2, double> & b,
   const tarch::la::Vector<2, double> & c,
   const tarch::la::Vector<2, double> & d,
   tarch::la::Vector<2, double>       & intersectionPoint );

}} // namespace precice, query

#endif /* PRECICE_QUERY_GEOMETRYCOMPUTATIONS_HPP_ */
