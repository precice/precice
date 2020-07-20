#include "FindClosestTriangle.hpp"
#include <Eigen/Dense>
// #include <Eigen/Core>
#include "math/barycenter.hpp"
#include "math/differences.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace query {

FindClosestTriangle::FindClosestTriangle(
    const Eigen::VectorXd &searchPoint)
    : _searchPoint(searchPoint),
      _vectorToProjectionPoint(Eigen::VectorXd::Constant(_searchPoint.size(), std::numeric_limits<double>::max())),
      _parametersProjectionPoint({_shortestDistance, _shortestDistance, _shortestDistance})
{
}

const Eigen::VectorXd &FindClosestTriangle::getSearchPoint() const
{
  return _searchPoint;
}

bool FindClosestTriangle::hasFound() const
{
  return _closestTriangle != nullptr;
}

double FindClosestTriangle::getEuclidianDistance() const
{
  return _shortestDistance;
}

mesh::Triangle &FindClosestTriangle::getClosestTriangle()
{
  PRECICE_ASSERT(_closestTriangle != nullptr);
  return *_closestTriangle;
}

const Eigen::VectorXd &FindClosestTriangle::getVectorToProjectionPoint() const
{
  return _vectorToProjectionPoint;
}

double FindClosestTriangle::getProjectionPointParameter(int index) const
{
  return _parametersProjectionPoint[index];
}

void FindClosestTriangle::find(
    mesh::Triangle &triangle)
{
  using Eigen::Vector2d;
  using Eigen::Vector3d;

  auto ret = math::barycenter::calcBarycentricCoordsForTriangle(
      triangle.vertex(0).getCoords(),
      triangle.vertex(1).getCoords(),
      triangle.vertex(2).getCoords(),
      triangle.getNormal(),
      _searchPoint);
  PRECICE_ASSERT(ret.barycentricCoords.size() == 3);

  // Determine from barycentric coordinates, if point is inside triangle
  const bool inside = not(ret.barycentricCoords.array() < -math::NUMERICAL_ZERO_DIFFERENCE).any();

  // If inside, compute distance to triangle and evtl. store distance
  if (inside) {
    Vector3d distanceVector = ret.projected - _searchPoint;
    double   distance       = distanceVector.norm();
    if (_shortestDistance > distance) {
      _shortestDistance             = distance;
      _vectorToProjectionPoint      = distanceVector;
      _parametersProjectionPoint[0] = ret.barycentricCoords(0);
      _parametersProjectionPoint[1] = ret.barycentricCoords(1);
      _parametersProjectionPoint[2] = ret.barycentricCoords(2);
      _closestTriangle              = &triangle;
    }
  }
}

} // namespace query
} // namespace precice
