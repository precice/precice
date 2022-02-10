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
    const Eigen::VectorXd &edgeA,
    const Eigen::VectorXd &edgeB,
    const Eigen::VectorXd &location)
{
  using Eigen::Vector2d;
  using Eigen::VectorXd;

  const int dimensions = edgeA.size();
  PRECICE_ASSERT(dimensions == edgeB.size(), "A and B need to have the same dimensions.", dimensions, edgeB.size());
  PRECICE_ASSERT(dimensions == location.size(), "A and the location need to have the same dimensions.", dimensions, location.size());
  PRECICE_ASSERT((dimensions == 2) || (dimensions == 3), dimensions);

  Vector2d barycentricCoords;
  VectorXd ab, ap;
  double   lenAb, lenProjected;

  // constant per edge
  ab    = edgeB - edgeA;
  lenAb = sqrt(ab.dot(ab));

  // varying per point
  ap           = location - edgeA;
  lenProjected = ap.dot(ab / lenAb);

  barycentricCoords(1) = lenProjected / lenAb;
  barycentricCoords(0) = 1 - barycentricCoords(1);

  return barycentricCoords;
}

Eigen::Vector3d calcBarycentricCoordsForTriangle(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &location)
{
  using Eigen::Vector3d;

  const int dimensions = a.size();
  PRECICE_ASSERT(dimensions == 3, dimensions);
  PRECICE_ASSERT(dimensions == b.size(), "A and B need to have the same dimensions.", dimensions, b.size());
  PRECICE_ASSERT(dimensions == c.size(), "A and C need to have the same dimensions.", dimensions, c.size());
  PRECICE_ASSERT(dimensions == location.size(), "A and the location need to have the same dimensions.", dimensions, location.size());

  Vector3d ab, ac, ap, n, barycentricCoords;
  double   scaleFactor;

  // constant per triangle
  ab          = b - a;
  ac          = c - a;
  n           = ab.cross(ac);
  scaleFactor = 1.0 / n.dot(n);

  // varying per point
  ap = location - a;

  barycentricCoords(2) = n.dot(ab.cross(ap)) * scaleFactor;
  barycentricCoords(1) = n.dot(ap.cross(ac)) * scaleFactor;
  barycentricCoords(0) = 1 - barycentricCoords(1) - barycentricCoords(2);

  return barycentricCoords;
}

} // namespace barycenter
} // namespace math
} // namespace precice
