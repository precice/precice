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

BarycentricCoordsAndProjected calcBarycentricCoordsForEdge(
    const Eigen::VectorXd &edgeA,
    const Eigen::VectorXd &edgeB,
    const Eigen::VectorXd &edgeNormal,
    const Eigen::VectorXd &location)
{
  using Eigen::Vector2d;
  using Eigen::Vector3d;
  using Eigen::VectorXd;

  const int dimensions = edgeA.size();
  PRECICE_ASSERT(dimensions == edgeB.size() && dimensions == edgeNormal.size() && dimensions == location.size(),
                 "The inputs need to have the same dimensions.");
  PRECICE_ASSERT((dimensions == 2) || (dimensions == 3), dimensions);

  Vector2d barycentricCoords;
  VectorXd projected = VectorXd::Zero(dimensions);
  Vector2d a, b, ab, c, d;
  bool     collinear = false;

  if (dimensions == 2) {
    // Get parameters for parametric edge representation: p(s) = a + s(b-a)
    a  = edgeA.head<2>();
    b  = edgeB.head<2>();
    ab = b;
    ab -= a;
    // Same for intersecting normal from searchpoint: q(t) = c + t(d - c)
    c = location.head<2>();
    d = location.head<2>();
    d += edgeNormal;
    collinear = math::geometry::collinear(a, b, c);
    if (collinear) {
      // From p(s) = a + s(b-a) we get: s = (p(s) - a) / (b-a)
      int iMax;
      ab.cwiseAbs().maxCoeff(&iMax);
      PRECICE_ASSERT(!math::equals(ab(iMax), 0.0));
      barycentricCoords[0] = (location(iMax) - a(iMax)) / ab(iMax);
      barycentricCoords[1] = 1.0 - barycentricCoords[0];
      projected            = location.head<2>();
      std::swap(barycentricCoords(0), barycentricCoords(1));
      return {barycentricCoords, projected};
    }
  } else { // 3D
    PRECICE_ASSERT(dimensions == 3, dimensions);
    // Get parameters for parametric triangle representation: p(s) = a + s(b-a)
    Vector3d a3D  = edgeA;
    Vector3d b3D  = edgeB;
    Vector3d c3D  = location;
    Vector3d ab3D = b3D - a3D;
    Vector3d ac3D = c3D - a3D;
    collinear     = math::geometry::collinear(a3D, b3D, c3D);
    if (collinear) {
      // From p(s) = a + s(b-a) we get: s = (p(s) - a) / (b-a)
      int iMax;
      ab3D.cwiseAbs().maxCoeff(&iMax);
      PRECICE_ASSERT(!math::equals(ab3D[iMax], 0.0));
      barycentricCoords[0] =
          (location(iMax) - edgeA(iMax)) / ab3D(iMax);
      barycentricCoords[1] = 1.0 - barycentricCoords[0];
      projected            = location;
      std::swap(barycentricCoords(0), barycentricCoords(1));
      return {barycentricCoords, projected};
    }
    // Project parameters to 2D, where the projection plane is determined from
    // the normal direction, in order to prevent "faulty" projections.
    Vector3d normal;
    normal = ab3D.cross(ac3D);
    int indexToRemove;
    normal.cwiseAbs().maxCoeff(&indexToRemove);
    int indices[2];
    if (indexToRemove == 0) {
      indices[0] = 1;
      indices[1] = 2;
    } else if (indexToRemove == 1) {
      indices[0] = 0;
      indices[1] = 2;
    } else {
      PRECICE_ASSERT(indexToRemove == 2, indexToRemove);
      indices[0] = 0;
      indices[1] = 1;
    }
    a << a3D[indices[0]], a3D[indices[1]];
    b << b3D[indices[0]], b3D[indices[1]];
    ab << ab3D[indices[0]], ab3D[indices[1]];
    c << c3D[indices[0]], c3D[indices[1]];
    // 3D normal might be projected out, hence, compute new 2D edge normal
    Vector2d normal2D(-1.0 * ab(1), ab(0));
    d = c + normal2D;
  }

  // Compute denominator for solving 2x2 equation system
  double D = a(0) * (d(1) - c(1)) + b(0) * (c(1) - d(1)) + d(0) * ab(1) - c(0) * ab(1);
  PRECICE_ASSERT(not math::equals(D, 0.0), a, b, c, d, ab); // D == 0 would imply "normal // edge"

  // Compute triangle segment parameter s, which is in [0, 1], if the
  // intersection point is within the triangle.
  barycentricCoords[0] = (a(0) * (d(1) - c(1)) +
                          c(0) * (a(1) - d(1)) +
                          d(0) * (c(1) - a(1))) /
                         D;
  barycentricCoords[1] = 1.0 - barycentricCoords[0];

  // Compute coordinates of projected point, which has to be done dimension
  // dependent again:
  if (dimensions == 2) {
    projected = ab;                    // = b - a
    projected *= barycentricCoords[0]; // = bary0 * (b - a)
    projected += a;                    // = a + bary0 * (b - a)
  } else {
    projected = edgeB;                 // = b
    projected -= edgeA;                // = b - a
    projected *= barycentricCoords[0]; // = bary0 * (b - a)
    projected += edgeA;                // = a + bary0 * (b - a)
  }

  std::swap(barycentricCoords(0), barycentricCoords(1));
  return {barycentricCoords, projected};
}

BarycentricCoordsAndProjected calcBarycentricCoordsForTriangle(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &normal,
    const Eigen::VectorXd &location)
{
  using Eigen::Vector2d;
  using Eigen::Vector3d;

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
