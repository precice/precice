#pragma once

#include <Eigen/Core>
#include <array>
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace math {

/// Provides computational geometry operations.
namespace geometry {

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
 * @param[in] a First point of line segment 1
 * @param[in] b Second point of line segement 1
 * @param[in] c First point of line segment 2
 * @param[in] d Second point of line segment 2
 *
 * @return True, if interseting. False, otherwise
 */
bool segmentsIntersect(
    const Eigen::Ref<const Eigen::Vector2d> &a,
    const Eigen::Ref<const Eigen::Vector2d> &b,
    const Eigen::Ref<const Eigen::Vector2d> &c,
    const Eigen::Ref<const Eigen::Vector2d> &d,
    bool                                     countTouchingAsIntersection);

/**
 * @brief Determines the intersection point of two lines, if one exists.
 *
 * Works only for Dim2.
 * @param[in] a First point on first line.
 * @param[in] b Second point on first line.
 * @param[in] c First point on second line.
 * @param[in] d Second point on second line.
 * @param[out] intersectionPoint Point of intersection.
 * @return True, if an intersection point exists and could be determined.
 */
bool lineIntersection(
    const Eigen::Ref<const Eigen::Vector2d> &a,
    const Eigen::Ref<const Eigen::Vector2d> &b,
    const Eigen::Ref<const Eigen::Vector2d> &c,
    const Eigen::Ref<const Eigen::Vector2d> &d,
    Eigen::Ref<Eigen::Vector2d> &            intersectionPoint);

/**
 * @brief Determines the intersection point of a segment with a plane in 3D.
 *
 * @param[out] intersectionPoint Contains coordinates of intersection point,
 *        if return value equals TOUCHING or INTERSECTION.
 *
 * @return Returns type of intersection. One of
 *         - NO_INTERSECTION
 *         - INTERSECTION
 *         - TOUCHING
 *         - CONTAINED
 */
ResultConstants segmentPlaneIntersection(
    const Eigen::Vector3d &pointOnPlane,
    const Eigen::Vector3d &planeNormal,
    const Eigen::Vector3d &firstPointSegment,
    const Eigen::Vector3d &secondPointSegment,
    Eigen::Vector3d &      intersectionPoint);

/**
 * @brief Determines, if a point lies on the line segment defined by a, b
 *
 * Works for 2D and 3D.
 *
 * @param[in] a       First point of line segment
 * @param[in] b       Second point of line segment
 * @param[in] toCheck Point to be checked for "betweenness"
 *
 * @return True, if toCheck lies between a and b. False, otherwise
 */
template <typename VECTORA_T, typename VECTORB_T, typename VECTORC_T>
bool between(
    const VECTORA_T &a,
    const VECTORB_T &b,
    const VECTORC_T &toCheck);

/**
 * @brief Determines, if three points are collinear (on one line)
 *
 * Works for Dim2 and Dim3.
 *
 * @param[in] a First point to check
 * @param[in] b Second point to check
 * @param[in] c Third point to check
 *
 * @return True, if collinear. False, otherwise
 */
template <typename Derived>
bool collinear(
    const Eigen::MatrixBase<Derived> &a,
    const Eigen::MatrixBase<Derived> &b,
    const Eigen::MatrixBase<Derived> &c);

/**
 * @brief Determines, if two lines are parallel to each other.
 *
 * Works for Dim2 and Dim3.
 *
 * @param[in] a First point on first line.
 * @param[in] b Second point on first line.
 * @param[in] c First point on second line.
 * @param[in] d Second point on second line.
 * @return True, if the two lines are parallel to each other.
 */
template <typename Derived>
static bool parallel(
    const Eigen::MatrixBase<Derived> &a,
    const Eigen::MatrixBase<Derived> &b,
    const Eigen::MatrixBase<Derived> &c,
    const Eigen::MatrixBase<Derived> &d);

/**
 * @brief Computes the signed area of a triangle in 2D.
 *
 * The area is negative, if the points a, b, c are given in
 * clockwise ordering, otherwise positive.
 */
double triangleArea(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c);

/// Computes the (unsigned) area of a triangle in 3D.
double tetraVolume(
    const Eigen::Vector3d &a,
    const Eigen::Vector3d &b,
    const Eigen::Vector3d &c,
    const Eigen::Vector3d &d);

/// Projects a 3D vector to a 2D one by removing one dimension.
Eigen::Vector2d projectVector(
    const Eigen::Vector3d &vector,
    const int              indexDimensionToRemove);

/**
 * @brief Tests if a vertex is contained in a triangle.
 *
 * @return One of:
 *         - CONTAINED
 *         - NOT_CONTAINED
 *         - TOUCHING
 */
int containedInTriangle(
    const Eigen::Vector2d &triangleVertex0,
    const Eigen::Vector2d &triangleVertex1,
    const Eigen::Vector2d &triangleVertex2,
    const Eigen::Vector2d &testPoint);

/**
 * @brief Tests, if a vertex is contained in a hyperrectangle.
 *
 * @return One of:
 *         - CONTAINED
 *         - NOT_CONTAINED
 *         - TOUCHING
 */
template <class Derived>
int containedInHyperrectangle(
    const Eigen::MatrixBase<Derived> &sidelengths,
    const Eigen::MatrixBase<Derived> &center,
    const Eigen::MatrixBase<Derived> &testPoint);

// --------------------------------------------------------- HEADER DEFINITIONS

template <typename VECTORA_T, typename VECTORB_T, typename VECTORC_T>
bool between(
    const VECTORA_T &a,
    const VECTORB_T &b,
    const VECTORC_T &toCheck)
{
  if (not collinear(a, b, toCheck)) {
    return false;
  }

  if (a(0) != b(0)) {
    return (math::greaterEquals(toCheck(0), a(0)) && math::greaterEquals(b(0), toCheck(0))) ||
           (math::greaterEquals(a(0), toCheck(0)) && math::greaterEquals(toCheck(0), b(0)));
  } else {
    return (math::greaterEquals(toCheck(1), a(1)) && math::greaterEquals(b(1), toCheck(1))) ||
           (math::greaterEquals(a(1), toCheck(1)) && math::greaterEquals(toCheck(1), b(1)));
  }
}

template <typename Derived>
bool collinear(
    const Eigen::MatrixBase<Derived> &a,
    const Eigen::MatrixBase<Derived> &b,
    const Eigen::MatrixBase<Derived> &c)
{
  PRECICE_ASSERT(a.size() == b.size(), a.size(), b.size());
  PRECICE_ASSERT(a.size() == c.size(), a.size(), c.size());
  double triangleOutline = (b - a).norm() + (c - b).norm() + (a - c).norm();
  if (math::equals(triangleArea(a, b, c) / triangleOutline, 0.0)) {
    return true;
  }
  return false;
}

template <class Derived>
bool parallel(
    const Eigen::MatrixBase<Derived> &a,
    const Eigen::MatrixBase<Derived> &b,
    const Eigen::MatrixBase<Derived> &c,
    const Eigen::MatrixBase<Derived> &d)
{
  if (math::equals(triangleArea(a, b, c), 0.0) and math::equals(triangleArea(a, b, d), 0.0))
    return true;

  return false;
}

template <class Derived>
int containedInHyperrectangle(
    const Eigen::MatrixBase<Derived> &sidelengths,
    const Eigen::MatrixBase<Derived> &center,
    const Eigen::MatrixBase<Derived> &testPoint)
{
  int             dim      = sidelengths.size();
  Eigen::VectorXd toCenter = testPoint - center;
  toCenter                 = toCenter.cwiseAbs();

  double diff     = 0.0;
  bool   touching = false;
  for (int i = 0; i < dim; i++) {
    diff = 0.5 * sidelengths(i) - toCenter(i);
    if (math::greater(0.0, diff)) {
      return NOT_CONTAINED;
    }
    if (math::equals(diff, 0.0)) {
      touching = true;
    }
  }
  if (touching) {
    return TOUCHING;
  }
  return CONTAINED;
}

struct ConvexityResult {
  bool               convex;
  std::array<int, 4> vertexOrder;
};

ConvexityResult isConvexQuad(std::array<Eigen::VectorXd, 4> coords);

} // namespace geometry
} // namespace math
} // namespace precice
