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
  double lenAb, lenProjected;

  // constant per edge
  ab = edgeB - edgeA;
  lenAb = sqrt(ab.dot(ab));

  // varying per point
  ap = location - edgeA;
  lenProjected = ap.dot(ab / lenAb);
  
  barycentricCoords(1) = lenProjected / lenAb;
  barycentricCoords(0) = 1 - barycentricCoords(1);

  return barycentricCoords;
}

Eigen::Vector3d calcBarycentricCoordsForTriangle(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &normal,
    const Eigen::VectorXd &location)
{
  using Eigen::Vector2d;
  using Eigen::Vector3d;

  const int dimensions = a.size();
  PRECICE_ASSERT(dimensions == 3, dimensions);
  PRECICE_ASSERT(dimensions == b.size(), "A and B need to have the same dimensions.", dimensions, b.size());
  PRECICE_ASSERT(dimensions == c.size(), "A and C need to have the same dimensions.", dimensions, c.size());
  PRECICE_ASSERT(dimensions == normal.size(), "A and the normal need to have the same dimensions.", dimensions, normal.size());
  PRECICE_ASSERT(dimensions == location.size(), "A and the location need to have the same dimensions.", dimensions, location.size());

  // Parametric representation for triangle plane:
  // (x, y, z) * normal = d
  const double d = normal.dot(a);

  // Parametric description of line from searchpoint orthogonal to triangle:
  // location + t * normal = x     (where t is parameter)
  // Determine t such that x lies on triangle plane:
  const double t = d - location.dot(normal) / normal.dot(normal);

  // Compute projected point with parameter t:
  Vector3d projected = normal;
  projected *= t;
  projected += location;

  // Project everything to 2D
  int iMax;
  normal.cwiseAbs().maxCoeff(&iMax);
  int indices[2];
  if (iMax == 0) {
    indices[0] = 1;
    indices[1] = 2;
  } else if (iMax == 1) {
    indices[0] = 0;
    indices[1] = 2;
  } else {
    PRECICE_ASSERT(iMax == 2, iMax);
    indices[0] = 0;
    indices[1] = 1;
  }
  Vector2d a2D(a[indices[0]],
               a[indices[1]]);
  Vector2d b2D(b[indices[0]],
               b[indices[1]]);
  Vector2d c2D(c[indices[0]],
               c[indices[1]]);
  Vector2d projected2D(projected[indices[0]], projected[indices[1]]);
  // Compute barycentric coordinates by solving linear 3x3 system
  Vector3d                    rhs(projected2D(0), projected2D(1), 1);
  Eigen::Matrix<double, 3, 3> A;
  A << a2D(0), b2D(0), c2D(0),
      a2D(1), b2D(1), c2D(1),
      1, 1, 1;
  Vector3d barycentricCoords = A.colPivHouseholderQr().solve(rhs);

  return {barycentricCoords, projected};
}

} // namespace barycenter
} // namespace math
} // namespace precice
