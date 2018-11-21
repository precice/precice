#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include "math/geometry.hpp"

namespace precice
{
namespace math
{
namespace barycenter
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
 *
 *  @param edgeA point A of the edge AB
 *  @param edgeNormal the normal of the edge
 *  @param location the location to compute the barycentric coordinates for
 *
 * @note Methodology of book "Computational Geometry", Joseph O' Rourke, Chapter 7.2
 */
BarycentricCoordsAndProjected calcBarycentricCoordsForEdge(
    const Eigen::VectorXd &edgeA,
    const Eigen::VectorXd &edgeB,
    const Eigen::VectorXd &edgeNormal,
    const Eigen::VectorXd &location);

/** Takes the corner vertices of a triangle and its norm.
 *  It then calculates the projection of a location vector and generates the barycentric coordinates for the corner points.
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

/** Takes the corner vertices of a quad and its norm.
 *  It then calculates the projection of a location vector and generates the barycentric coordinates for the corner points.
 *
 *  @param a point A of the quad ABCD
 *  @param b point B of the quad ABCD
 *  @param c point C of the quad ABCD
 *  @param d point D of the quad ABCD
 *  @param normal the normal of the quad
 *  @param location the location to compute the barycentric coordinates for
 *
 *   @todo: Interpolation on quads is currently not implemented
 */
BarycentricCoordsAndProjected calcBarycentricCoordsForQuad(
    const Eigen::VectorXd &a,
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &c,
    const Eigen::VectorXd &d,
    const Eigen::VectorXd &normal,
    const Eigen::VectorXd &location);

} // namespace barycenter
} // namespace math
} // namespace precice
