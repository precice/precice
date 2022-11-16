#include "math/barycenter.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <utility>
#include "math/differences.hpp"
#include "math/geometry.hpp"
#include "utils/assertion.hpp"

namespace precice::math::barycenter {

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

  // Let AB be the edge and U the input point. We compute the projection P of U on the edge.
  // To find P, start from A and move from dot(AU, AB) / |AB| along the AB direction.
  // Divide again by |AB| to find the barycentric coordinate of P relative to U
  // This means we just need to compute dot(AU, AB) / dot(AB, AB)
  ab = b - a;
  au = u - a;

  barycentricCoords(1) = au.dot(ab) / ab.dot(ab);
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

    ab         = b - a;
    ac         = c - a;
    n          = ab.cross(ac);
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
    scaleFactor = 1.0 / twiceArea;

    barycentricCoords(0) = crossProduct2D(ub, uc) * scaleFactor;
    barycentricCoords(1) = crossProduct2D(uc, ua) * scaleFactor;
    barycentricCoords(2) = 1 - barycentricCoords(0) - barycentricCoords(1);
  }

  return barycentricCoords;
}

Eigen::Vector4d calcBarycentricCoordsForTetrahedron(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &d,
    const Eigen::VectorXd &u)
{
  using Eigen::Vector3d;
  using Eigen::Vector4d;

  const int dimensions = a.size();

  PRECICE_ASSERT(dimensions == 3, dimensions);
  PRECICE_ASSERT(dimensions == b.size(), "A and B need to have the same dimensions.", dimensions, b.size());
  PRECICE_ASSERT(dimensions == c.size(), "A and C need to have the same dimensions.", dimensions, c.size());
  PRECICE_ASSERT(dimensions == d.size(), "A and D need to have the same dimensions.", dimensions, d.size());
  PRECICE_ASSERT(dimensions == u.size(), "A and the point need to have the same dimensions.", dimensions, u.size());

  Vector4d barycentricCoords;

  // Varying per point
  Vector3d au = u - a;
  Vector3d du = u - d;
  Vector3d cu = u - c;

  // Necessary to compute triangles
  Vector3d ab = b - a;
  Vector3d ac = c - a;
  Vector3d ad = d - a;
  Vector3d bc = c - b;

  // Triangles
  Vector3d abc = ab.cross(bc);
  Vector3d abd = ab.cross(-ad);
  Vector3d acd = ac.cross(ad);

  auto volume = abc.dot(ad);

  barycentricCoords(3) = abc.dot(au) / volume;
  barycentricCoords(2) = abd.dot(du) / volume;
  barycentricCoords(1) = acd.dot(cu) / volume;
  barycentricCoords(0) = 1 - barycentricCoords(3) - barycentricCoords(2) - barycentricCoords(1);

  return barycentricCoords;
}

} // namespace precice::math::barycenter
