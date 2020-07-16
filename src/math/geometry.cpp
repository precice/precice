#include "geometry.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math.h>
#include <stdlib.h>
#include "math/math.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace math {
namespace geometry {

bool segmentsIntersect(
    const Eigen::Ref<const Eigen::Vector2d> &a,
    const Eigen::Ref<const Eigen::Vector2d> &b,
    const Eigen::Ref<const Eigen::Vector2d> &c,
    const Eigen::Ref<const Eigen::Vector2d> &d,
    bool                                     countTouchingAsIntersection)
{
  PRECICE_ASSERT(a.size() == 2, a.size());
  PRECICE_ASSERT(b.size() == 2, b.size());
  PRECICE_ASSERT(c.size() == 2, c.size());
  PRECICE_ASSERT(d.size() == 2, d.size());

  if (countTouchingAsIntersection) {
    if (between(a, b, c)) {
      return true;
    } else if (between(a, b, d)) {
      return true;
    } else if (between(c, d, a)) {
      return true;
    } else if (between(c, d, b)) {
      return true;
    }
  }

  double abc = triangleArea(a, b, c);
  double abd = triangleArea(a, b, d);
  double cda = triangleArea(c, d, a);
  double cdb = triangleArea(c, d, b);

  double circABC = (a - b).norm() + (b - c).norm() + (c - a).norm();
  double circABD = (a - b).norm() + (b - d).norm() + (d - a).norm();
  double circCDA = (c - d).norm() + (d - a).norm() + (a - c).norm();
  double circCDB = (c - d).norm() + (d - b).norm() + (b - c).norm();

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
  using utils::xOR;

  if (xOR(std::abs(abc) <= math::NUMERICAL_ZERO_DIFFERENCE,
          std::abs(abd) <= math::NUMERICAL_ZERO_DIFFERENCE)) {
    return false;
  }
  if (xOR(std::abs(cda) <= math::NUMERICAL_ZERO_DIFFERENCE,
          std::abs(cdb) <= math::NUMERICAL_ZERO_DIFFERENCE)) {
    return false;
  }

  // Check, whether the segments are intersecting in the real sense.
  bool isFirstSegmentBetween  = xOR(abc > -math::NUMERICAL_ZERO_DIFFERENCE,
                                   abd > -math::NUMERICAL_ZERO_DIFFERENCE);
  bool isSecondSegmentBetween = xOR(cda > -math::NUMERICAL_ZERO_DIFFERENCE,
                                    cdb > -math::NUMERICAL_ZERO_DIFFERENCE);

  return isFirstSegmentBetween && isSecondSegmentBetween;
}

bool lineIntersection(
    const Eigen::Ref<const Eigen::Vector2d> &a,
    const Eigen::Ref<const Eigen::Vector2d> &b,
    const Eigen::Ref<const Eigen::Vector2d> &c,
    const Eigen::Ref<const Eigen::Vector2d> &d,
    Eigen::Ref<Eigen::Vector2d> &            intersectionPoint)
{
  // Compute denominator for solving 2x2 equation system
  double D = a(0) * (d(1) - c(1)) +
             b(0) * (c(1) - d(1)) +
             d(0) * (b(1) - a(1)) -
             c(0) * (a(1) - b(1));

  // If D==0, the two lines are parallel
  if (math::equals(D, 0.0)) {
    return false;
  }

  // Compute line segment parameter s, which is in [0, 1], if the
  // intersection point is within the segment.
  double s = (a(0) * (d(1) - c(1)) + c(0) * (a(1) - d(1)) + d(0) * (c(1) - a(1))) / D;

  intersectionPoint = a + s * (b - a);
  return true;
}

ResultConstants segmentPlaneIntersection(
    const Eigen::Vector3d &pointOnPlane,
    const Eigen::Vector3d &planeNormal,
    const Eigen::Vector3d &firstPointSegment,
    const Eigen::Vector3d &secondPointSegment,
    Eigen::Vector3d &      intersectionPoint)
{
  // Methodology of "Computation Geometry in C", Joseph O'Rourke, Chapter 7.3.1

  // The plane is represented by: p(x,y,z) dot planeNormal = d  (1)
  // We need to compute d by using pointOnPlane as p(x,y,z):
  double d = pointOnPlane.dot(planeNormal);

  // The segment is represented by:
  // firstPointSegment + (secondPointSegment - firstPointSegement) * t (2)
  // We insert (2) into (1) and reorder for t. The resulting equations reads
  // t = (d - firstPointSegment dot planeNormal) /
  //     ((secondPointSegment - firstPointSegment) dot planeNormal). (3)
  //
  // If the denominator of (3) equals 0, the segment is parallel to the plane.
  // If, in addition, the denominator of (3) equals 0, the segment lies in the
  // plane. Otherwise we can compute a point of intersection.
  Eigen::Vector3d segmentVec(secondPointSegment - firstPointSegment);
  double          nominator   = d - firstPointSegment.dot(planeNormal);
  double          denominator = segmentVec.dot(planeNormal);
  if (math::equals(denominator, 0.0)) {
    if (math::equals(nominator, 0.0)) {
      return CONTAINED;
    } else {
      return NO_INTERSECTION;
    }
  }
  double t = nominator / denominator;

  // If t is larger than 1 or smaller than zero, the intersection is not within
  // the line segment
  if (math::greater(t, 1.0) || math::greater(0.0, t)) {
    return NO_INTERSECTION;
  }

  intersectionPoint = firstPointSegment + segmentVec * t;

  // If t equals 1 or 0, the segment is just touching the plane, otherwise, a
  // real intersection is happening.
  if (math::equals(t, 0.0) || math::equals(t, 1.0)) {
    return TOUCHING;
  }

  return INTERSECTION;
}

double triangleArea(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c)
{
  PRECICE_ASSERT(a.size() == b.size(), a.size(), b.size());
  PRECICE_ASSERT(b.size() == c.size(), b.size(), c.size());
  if (a.size() == 2) {
    Eigen::Vector2d A = b;
    A -= a;
    Eigen::Vector2d B = c;
    B -= a;
    return 0.5 * (A(0) * B(1) - A(1) * B(0));
  } else {
    PRECICE_ASSERT(a.size() == 3, a.size());
    Eigen::Vector3d A = b;
    A -= a;
    Eigen::Vector3d B = c;
    B -= a;
    return 0.5 * A.cross(B).norm();
  }
}

double tetraVolume(
    const Eigen::Vector3d &a,
    const Eigen::Vector3d &b,
    const Eigen::Vector3d &c,
    const Eigen::Vector3d &d)
{
  return std::abs((a - d).dot((b - d).cross(c - d))) / 6.0;
}

Eigen::Vector2d projectVector(
    const Eigen::Vector3d &vector,
    int                    indexDimensionToRemove)
{
  PRECICE_ASSERT(indexDimensionToRemove >= 0);
  PRECICE_ASSERT(indexDimensionToRemove < 3);
  Eigen::Vector2d projectedVector;
  int             reducedDim = 0;
  for (int dim = 0; dim < 3; dim++) {
    if (indexDimensionToRemove == dim) {
      continue;
    }
    projectedVector(reducedDim) = vector(dim);
    reducedDim++;
  }
  return projectedVector;
}

int containedInTriangle(
    const Eigen::Vector2d &triangleVertex0,
    const Eigen::Vector2d &triangleVertex1,
    const Eigen::Vector2d &triangleVertex2,
    const Eigen::Vector2d &testPoint)
{
  double area0 = triangleArea(triangleVertex0, triangleVertex1, testPoint);
  double area1 = triangleArea(triangleVertex1, triangleVertex2, testPoint);
  double area2 = triangleArea(triangleVertex2, triangleVertex0, testPoint);

  double length01 = (triangleVertex0 - triangleVertex1).norm();
  double length12 = (triangleVertex1 - triangleVertex2).norm();
  double length20 = (triangleVertex2 - triangleVertex0).norm();
  double length1t = (triangleVertex1 - testPoint).norm();
  double length2t = (triangleVertex2 - testPoint).norm();
  double length0t = (triangleVertex0 - testPoint).norm();
  double scale0   = length01 + length1t + length0t;
  double scale1   = length12 + length2t + length1t;
  double scale2   = length20 + length0t + length2t;
  area0 /= scale0;
  area1 /= scale1;
  area2 /= scale2;

  int sign0       = math::sign(area0);
  int sign1       = math::sign(area1);
  int sign2       = math::sign(area2);
  int absSumSigns = std::abs(sign0 + sign1 + sign2);
  if (absSumSigns == 3) {
    // In triangle
    return CONTAINED;
  } else if (absSumSigns == 2) {
    // On triangle edge
    return TOUCHING;
  } else if (std::abs(sign0) + std::abs(sign1) + std::abs(sign2) == 1) {
    // On triangle vertex
    return TOUCHING;
  }
  // Outside of triangle
  return NOT_CONTAINED;
}

} // namespace geometry
} // namespace math
} // namespace precice
