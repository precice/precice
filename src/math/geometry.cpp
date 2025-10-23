#include "geometry.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "logging/Logger.hpp"
#include "math/math.hpp"
#include "utils/Helpers.hpp"
#include "utils/assertion.hpp"

namespace precice::math::geometry {

logging::Logger _log("math::geometry");

bool lineIntersection(
    const Eigen::Ref<const Eigen::Vector2d> &a,
    const Eigen::Ref<const Eigen::Vector2d> &b,
    const Eigen::Ref<const Eigen::Vector2d> &c,
    const Eigen::Ref<const Eigen::Vector2d> &d,
    Eigen::Ref<Eigen::Vector2d>             &intersectionPoint)
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
    Eigen::Vector3d       &intersectionPoint)
{
  // Methodology of "Computation Geometry in C", Joseph O'Rourke, Chapter 7.3.1

  // The plane is represented by: p(x,y,z) dot planeNormal = d  (1)
  // We need to compute d by using pointOnPlane as p(x,y,z):
  double d = pointOnPlane.dot(planeNormal);

  // The segment is represented by:
  // firstPointSegment + (secondPointSegment - firstPointSegment) * t (2)
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
    return 0.5 * std::fabs(A(0) * B(1) - A(1) * B(0));
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

ConvexityResult isConvexQuad(std::array<Eigen::VectorXd, 4> coords)
{
  /*
    All points need to be projected into a new plane with only 2 coordinates, x' and y'. These are used to check
    the convexity of the quad. These new coordinates are stored in 'coords'.
  */
  PRECICE_ASSERT(std::all_of(coords.cbegin(), coords.cend(),
                             [](const auto &v) { return v.size() == 2 || v.size() == 3; }),
                 "This only works in 2D or 3D.");

  // Check that quad is planar and project to a new plane
  if (coords[0].size() == 3) {

    Eigen::Vector3d coordOrigin; // Origin point for the transformation of points onto the new plane
    coordOrigin = coords[0];

    // Normal of the plane of first three points in the list of vertices
    Eigen::Vector3d e_1          = coords[1] - coordOrigin;
    Eigen::Vector3d e_2          = coords[2] - coordOrigin;
    Eigen::Vector3d normalVector = e_1.cross(e_2);

    PRECICE_CHECK(math::equals(normalVector.dot(coords[3] - coordOrigin), 0.0),
                  "Non-planar quads are not supported. The vertex coordinates are: {}.", coords);

    // Transform Coordinates - coord[0] is the origin
    for (int i = 0; i < 4; i++) {
      Eigen::Vector3d coordinateDifference = coords[i] - coordOrigin;
      coords[i][0]                         = e_1.dot(coordinateDifference);
      coords[i][1]                         = e_2.dot(coordinateDifference);
      coords[i][2]                         = normalVector.dot(coordinateDifference);
    }
  }
  if (coords[0].size() == 2) {
    Eigen::Vector2d coordOrigin; // Origin point for the transformation of points onto the new plane
    coordOrigin = coords[0];

    Eigen::Vector2d e_1 = coords[1] - coordOrigin;
    Eigen::Vector2d e_2 = coords[2] - coordOrigin;

    // Transform Coordinates - coord[0] is the origin
    for (int i = 0; i < 4; i++) {
      Eigen::Vector2d coordinateDifference = coords[i] - coordOrigin;
      coords[i][0]                         = e_1.dot(coordinateDifference);
      coords[i][1]                         = e_2.dot(coordinateDifference);
    }
  }

  /*
  For the convex hull algorithm, the most left hand point regarding the x coordinate is chosen as the starting point.
  The algorithm moves in an anti-clockwise position, finding the most right hand coordinate from the
  previous most right hand point. The algorithm must find 3 other points in order to be a valid quad.
  */

  // First find point with smallest x coord. This point must be in the convex set then and is the starting point of gift wrapping algorithm
  int idLowestPoint = 0;
  for (int i = 1; i < 4; i++) {
    if (coords[i][0] < coords[idLowestPoint][0]) {
      idLowestPoint = i;
    }
  }

  // Found starting point. Add this as the first vertex in the convex hull.
  // current is the origin point => hull[0]
  int validVertexIDCounter = 0;             // Counts number of times a valid vertex is found. Must be 4 for a valid quad.
  int currentVertex        = idLowestPoint; // current valid vertex

  ConvexityResult result{};
  do {
    // Add current point to result
    result.vertexOrder[validVertexIDCounter] = currentVertex;

    // Next potential valid vertex
    int nextVertex = (currentVertex + 1) % 4; // remainder resets loop through vector of points
    for (int i = 0; i < 4; i++) {
      double y1  = coords[currentVertex][1] - coords[nextVertex][1];
      double y2  = coords[currentVertex][1] - coords[i][1];
      double x1  = coords[currentVertex][0] - coords[nextVertex][0];
      double x2  = coords[currentVertex][0] - coords[i][0];
      double val = y2 * x1 - y1 * x2;

      if (val > 0) {
        nextVertex = i;
      }
    }
    // Now nextVertex is the most anti-clockwise with respect to current
    // Set current as nextVertex for next iteration, so that nextVertex is added to
    // result 'hull'
    currentVertex = nextVertex;
    validVertexIDCounter++;
  } while (currentVertex != idLowestPoint); // While we don't come to first point

  // Ordering of quad is hull 0-1-2-3-0
  result.convex = (validVertexIDCounter == 4);

  return result;
}
} // namespace precice::math::geometry
