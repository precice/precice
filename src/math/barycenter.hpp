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
 *  @param p the point to compute the barycentric coordinates for
 *
 * @note Simple scalar projection approach, projected point is not actually calculated.
 */
Eigen::Vector2d calcBarycentricCoordsForEdge(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &p);

/** Takes the corner vertices of a triangle and a point in 3D space.
 *  Returns the barycentric coordinates for that point's projection onto the given triangle.
 *
 *  @param a point A of the triangle ABC
 *  @param b point B of the triangle ABC
 *  @param c point C of the triangle ABC
 *  @param normal the normal of the triangle
 *  @param location the location to compute the barycentric coordinates for
 *
 * @note
 * Methodology of book "Computational Geometry", Joseph O' Rourke, Chapter 7.3
 * with the barycentric coordinates method and real projection into 2D, instead
 * of outprojecting one coordinate
 */
BarycentricCoordsAndProjected calcBarycentricCoordsForTriangle(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &normal,
    const Eigen::VectorXd &location);

} // namespace barycenter
} // namespace math
} // namespace precice
