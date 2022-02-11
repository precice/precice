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

Eigen::Vector3d calcBarycentricCoordsForTriangle(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &u)
{
  using Eigen::Vector3d;

  const int dimensions = a.size();
  PRECICE_ASSERT(dimensions == 3, dimensions);
  PRECICE_ASSERT(dimensions == b.size(), "A and B need to have the same dimensions.", dimensions, b.size());
  PRECICE_ASSERT(dimensions == c.size(), "A and C need to have the same dimensions.", dimensions, c.size());
  PRECICE_ASSERT(dimensions == u.size(), "A and the point need to have the same dimensions.", dimensions, u.size());

  Vector3d ab, ac, au, n, barycentricCoords;
  double   scaleFactor;

  // constant per triangle
  ab          = b - a;
  ac          = c - a;
  n           = ab.cross(ac);
  scaleFactor = 1.0 / n.dot(n);

  // varying per point
  au = u - a;

  barycentricCoords(2) = n.dot(ab.cross(au)) * scaleFactor;
  barycentricCoords(1) = n.dot(au.cross(ac)) * scaleFactor;
  barycentricCoords(0) = 1 - barycentricCoords(1) - barycentricCoords(2);

  return barycentricCoords;
}

} // namespace barycenter
} // namespace math
} // namespace precice
