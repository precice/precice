#include "math/barycenter.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <utility>
#include "math/differences.hpp"
#include "math/geometry.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace math {
namespace barycenter {

Eigen::Vector2d calcBarycentricCoordsForEdge(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &u)
{
  using Eigen::Vector2d;
  using Eigen::VectorXd;

  const int dimensions = a.size();
  PRECICE_ASSERT(dimensions == b.size(), "A and B need to have the same dimensions.", dimensions, b.size());
  PRECICE_ASSERT(dimensions == u.size(), "A and the point need to have the same dimensions.", dimensions, u.size());
  PRECICE_ASSERT((dimensions == 2) || (dimensions == 3), dimensions);

  Vector2d barycentricCoords;
  VectorXd ab, au;
  double   lenAb, lenProjected;

  // constant per edge
  ab    = b - a;
  lenAb = sqrt(ab.dot(ab));

  // varying per point
  au           = u - a;
  lenProjected = au.dot(ab / lenAb);

  barycentricCoords(1) = lenProjected / lenAb;
  barycentricCoords(0) = 1 - barycentricCoords(1);

  return barycentricCoords;
}

static double crossProduct2D(const Eigen::Vector2d &u, const Eigen::Vector2d &v)
{
  return u(0) * v(1) - u(1) * v(0);
}

Eigen::Vector3d calcBarycentricCoordsForTriangle(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &u)
{
  using Eigen::Vector2d;
  using Eigen::Vector3d;

  const int dimensions = a.size();
  PRECICE_ASSERT(dimensions == b.size(), "A and B need to have the same dimensions.", dimensions, b.size());
  PRECICE_ASSERT(dimensions == c.size(), "A and C need to have the same dimensions.", dimensions, c.size());
  PRECICE_ASSERT(dimensions == u.size(), "A and the point need to have the same dimensions.", dimensions, u.size());

  Vector3d barycentricCoords;
  double   scaleFactor;

  // constant per triangle
  if (dimensions == 3) {
    Vector3d ab, ac, au, n;

    ab          = b - a;
    ac          = c - a;
    n           = ab.cross(ac);
    auto nDotN = n.dot(n);
    PRECICE_ASSERT(nDotN != 0, "It seems a degenerate triangle was sent.");
    scaleFactor = 1.0 / nDotN;

    // varying per point
    au = u - a;

    barycentricCoords(2) = n.dot(ab.cross(au)) * scaleFactor;
    barycentricCoords(1) = n.dot(au.cross(ac)) * scaleFactor;
    barycentricCoords(0) = 1 - barycentricCoords(1) - barycentricCoords(2);

  } else {
    // Dimension == 2 => No cross product
    Vector2d ab, ac, ub, uc, ua;
    ab = b - a;
    ac = c - a;
    ub = b - u;
    uc = c - u;
    ua = a - u;

    auto twiceArea = crossProduct2D(ab, ac);
    PRECICE_ASSERT(twiceArea != 0, "It seems a degenerate triangle was sent.");
    scaleFactor          = 1.0 / twiceArea;
    barycentricCoords(0) = crossProduct2D(ub, uc) * scaleFactor;
    barycentricCoords(1) = crossProduct2D(uc, ua) * scaleFactor;
    barycentricCoords(2) = 1 - barycentricCoords(0) - barycentricCoords(1);
  }

  return barycentricCoords;
}

} // namespace barycenter
} // namespace math
} // namespace precice
