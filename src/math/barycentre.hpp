#pragma once

#include <Eigen/Core>
#include "math/geometry.hpp"

namespace precice
{
namespace math
{
namespace barycentre
{

/// The result of calculating the barycentric coordinates.
struct BarycentricCoordsAndProjected {
  /// A vector of the n coefficients for n vertices
  Eigen::VectorXd barycentricCoords;
  /// The projected location vertex
  Eigen::VectorXd projected;
};

/** Takes the corner vertices of an edge and its norm.
 *  It then calculates the projection of a location vector and generates the barycentric coordinates for the corner points.
 */
template <class DerivedA, class DerivedB, class DerivedNorm, class DerivedLoc>
BarycentricCoordsAndProjected calcBarycentricCoordsForEdge(
        const Eigen::MatrixBase<DerivedA> &edgeA,
        const Eigen::MatrixBase<DerivedB> &edgeB,
        const Eigen::MatrixBase<DerivedNorm> &edgeNormal,
        const Eigen::MatrixBase<DerivedLoc> &location)
{
  using Eigen::Vector2d;
  using Eigen::Vector3d;
  Vector2d barycentricCoords;
  int      dimensions = edgeA.size();
  assertion(edgeA.size() == edgeB.size() == location.size());
  assertion((dimensions == 2) || (dimensions == 3), dimensions);
  Eigen::VectorXd projected = Eigen::VectorXd::Zero(dimensions);
  Vector2d        a, b, ab, c, d;
  bool            collinear = false;

  if (dimensions == 2) {
    // Get parameters for parametric edge representation: p(s) = a + s(b-a)
    a  = edgeA;
    b  = edgeB;
    ab = b;
    ab -= a;
    // Same for intersecting normal from searchpoint: q(t) = c + t(d - c)
    c = location;
    d = location;
    d += edgeNormal;
    collinear = math::geometry::collinear(a, b, c);
    if (collinear) {
      // From p(s) = a + s(b-a) we get: s = (p(s) - a) / (b-a)
      int iMax;
      ab.cwiseAbs().maxCoeff(&iMax);
      assertion(!math::equals(ab(iMax), 0.0));
      barycentricCoords[0] = (location(iMax) - a(iMax)) / ab(iMax);
      barycentricCoords[1] = 1.0 - barycentricCoords[0];
      projected            = location;
    }
  } else { // 3D
    assertion(dimensions == 3, dimensions);
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
      assertion(!math::equals(ab3D[iMax], 0.0));
      barycentricCoords[0] =
          (location(iMax) - edgeA(iMax)) / ab3D(iMax);
      barycentricCoords[1] = 1.0 - barycentricCoords[0];
      projected            = location;
    } else {
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
        assertion(indexToRemove == 2, indexToRemove);
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
  }

  if (not collinear) {
    // Compute denominator for solving 2x2 equation system
    double D = a(0) * (d(1) - c(1)) + b(0) * (c(1) - d(1)) + d(0) * ab(1) - c(0) * ab(1);
    assertion(not math::equals(D, 0.0), a, b, c, d, ab); // D == 0 would imply "normal // edge"

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
  }
  return {barycentricCoords, projected};
}

/** Takes the corner vertices of a triangle and its norm.
 *  It then calculates the projection of a location vector and generates the barycentric coordinates for the corner points.
 */
template <class DerivedA, class DerivedB, class DerivedC, class DerivedNorm, class DerivedLoc>
BarycentricCoordsAndProjected calcBarycentricCoordsForTriangle(
    const Eigen::MatrixBase<DerivedA> &   a,
    const Eigen::MatrixBase<DerivedB> &   b,
    const Eigen::MatrixBase<DerivedC> &   c,
    const Eigen::MatrixBase<DerivedNorm> &normal,
    const Eigen::MatrixBase<DerivedLoc> & location);

/** Takes the corner vertices of a quad and its norm.
 *  It then calculates the projection of a location vector and generates the barycentric coordinates for the corner points.
 */
template <class DerivedA, class DerivedB, class DerivedC, class DerivedD, class DerivedNorm, class DerivedLoc>
BarycentricCoordsAndProjected calcBarycentricCoordsForQuad(
    const Eigen::MatrixBase<DerivedA> &   a,
    const Eigen::MatrixBase<DerivedB> &   b,
    const Eigen::MatrixBase<DerivedC> &   c,
    const Eigen::MatrixBase<DerivedD> &   d,
    const Eigen::MatrixBase<DerivedNorm> &normal,
    const Eigen::MatrixBase<DerivedLoc> & location);

} // namespace barycentre
} // namespace math
} // namespace precice
