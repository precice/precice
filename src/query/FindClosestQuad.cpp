#include "FindClosestQuad.hpp"
#include <Eigen/Core>
#include "logging/LogMacros.hpp"
#include "math/barycenter.hpp"
#include "math/differences.hpp"
#include "mesh/Quad.hpp"
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace query {

FindClosestQuad::FindClosestQuad(
    const Eigen::VectorXd &searchPoint)
    : _searchPoint(searchPoint),
      _vectorToProjectionPoint(Eigen::VectorXd::Constant(_searchPoint.size(), std::numeric_limits<double>::max())),
      _parametersProjectionPoint({_shortestDistance, _shortestDistance, _shortestDistance, _shortestDistance})
{
}

const Eigen::VectorXd &FindClosestQuad::getSearchPoint() const
{
  return _searchPoint;
}

bool FindClosestQuad::hasFound() const
{
  return _closestQuad != nullptr;
}

double FindClosestQuad::getEuclidianDistance()
{
  return _shortestDistance;
}

mesh::Quad &FindClosestQuad::getClosestQuad()
{
  PRECICE_ASSERT(_closestQuad != nullptr);
  return *_closestQuad;
}

const Eigen::VectorXd &FindClosestQuad::getVectorToProjectionPoint() const
{
  return _vectorToProjectionPoint;
}

double FindClosestQuad::getProjectionPointParameter(
    int index) const
{
  return _parametersProjectionPoint[index];
}

void FindClosestQuad::find(
    mesh::Quad &quad)
{
  PRECICE_TRACE(quad.vertex(0).getCoords(), quad.vertex(1).getCoords(), quad.vertex(2).getCoords(), quad.vertex(3).getCoords());

  // Methodology of book "Computational Geometry", Joseph O' Rourke, Chapter 7.2
  auto &a    = quad.vertex(0).getCoords();
  auto &b    = quad.vertex(1).getCoords();
  auto &c    = quad.vertex(2).getCoords();
  auto &d    = quad.vertex(3).getCoords();
  auto &norm = quad.getNormal();

  auto ret = math::barycenter::calcBarycentricCoordsForQuad(a, b, c, d, norm, _searchPoint);
  PRECICE_ASSERT(ret.barycentricCoords.size() == 4);

  const bool inside = not(ret.barycentricCoords.array() < -math::NUMERICAL_ZERO_DIFFERENCE).any();

  // if valid, compute distance to triangle and evtl. store distance
  if (inside) {
    Eigen::VectorXd distanceVector = ret.projected;
    distanceVector -= _searchPoint;
    double distance = distanceVector.norm();
    if (_shortestDistance > distance) {
      _shortestDistance             = distance;
      _vectorToProjectionPoint      = distanceVector;
      _parametersProjectionPoint[0] = ret.barycentricCoords[0];
      _parametersProjectionPoint[1] = ret.barycentricCoords[1];
      _parametersProjectionPoint[2] = ret.barycentricCoords[2];
      _parametersProjectionPoint[3] = ret.barycentricCoords[3];
      _closestQuad                  = &quad;
    }
  }
}

} // namespace query
} // namespace precice
