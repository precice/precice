#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include "math/geometry.hpp"

/// Provides operations to calculate barycentric coordinates for a point's projection onto a primitive.
namespace precice::math::barycenter {

/** Takes the end vertices of an edge and a point in 2D or 3D space.
 *  Returns the barycentric coordinates for that point's projection onto the given edge.
 *
 *  @param a point A of the edge AB
 *  @param b point B of the edge AB
 *  @param u the point to compute the barycentric coordinates for
 *
 * @note Simple scalar projection approach, projected point in Cartesian coordinates is not actually calculated.
 */
Eigen::Vector2d calcBarycentricCoordsForEdge(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &u);

/** Takes the corner vertices of a triangle and a point in 3D space.
 *  Returns the barycentric coordinates for that point's projection onto the given triangle.
 *
 *  @param a point A of the triangle ABC
 *  @param b point B of the triangle ABC
 *  @param c point C of the triangle ABC
 *  @param u the point to compute the barycentric coordinates for
 *
 * @note This implements the efficient one-step algorithm (no separate projection) presented in
 *  _Computing the barycentric coordinates of a projected point_ by W. Heidrich (2005)
 *
 */
Eigen::Vector3d calcBarycentricCoordsForTriangle(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &u);

/** Takes the corner vertices of a tetrahedron and a point in 3D space.
 *  Returns the barycentric coordinates for that point's projection onto the given tetrahedron.
 *
 *  @param a point A of the tetrahedron ABCD
 *  @param b point B of the tetrahedron ABCD
 *  @param c point C of the tetrahedron ABCD
*  @param d point D of the tetrahedron ABCD
 *  @param u the point to compute the barycentric coordinates for
 *
 * @note This implements an efficient one-step algorithm (no separate projection)
 * described in Boris Martin's Master's thesis
 *
 */
Eigen::Vector4d calcBarycentricCoordsForTetrahedron(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &d,
    const Eigen::VectorXd &u);

} // namespace precice::math::barycenter
