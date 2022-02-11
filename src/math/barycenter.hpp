#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include "math/geometry.hpp"

namespace precice {
namespace math {
/// Provides operations to calculate barycentric coordinates for a point's projection onto a primitive.
namespace barycenter {

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

} // namespace barycenter
} // namespace math
} // namespace precice
