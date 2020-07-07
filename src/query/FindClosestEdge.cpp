#include "FindClosestEdge.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <array>
#include "logging/LogMacros.hpp"
#include "math/barycenter.hpp"
#include "math/differences.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace query {

FindClosestEdge::FindClosestEdge(
    const Eigen::VectorXd &searchPoint)
    : _searchPoint(searchPoint),
      _vectorToProjectionPoint(Eigen::VectorXd::Constant(searchPoint.size(), std::numeric_limits<double>::max())),
      _parametersProjectionPoint({_shortestDistance, _shortestDistance})
{
  PRECICE_ASSERT((_searchPoint.size() == 2) || (_searchPoint.size() == 3),
                 _searchPoint.size());
}

const Eigen::VectorXd &FindClosestEdge::getSearchPoint() const
{
  return _searchPoint;
}

bool FindClosestEdge::hasFound() const
{
  return _closestEdge != nullptr;
}

double FindClosestEdge::getEuclidianDistance()
{
  return _shortestDistance;
}

mesh::Edge &FindClosestEdge::getClosestEdge()
{
  PRECICE_ASSERT(_closestEdge != nullptr);
  return *_closestEdge;
}

const Eigen::VectorXd &FindClosestEdge::getVectorToProjectionPoint() const
{
  return _vectorToProjectionPoint;
}

double FindClosestEdge::getProjectionPointParameter(
    int index) const
{
  return _parametersProjectionPoint[index];
}

//double FindClosestEdge:: getFirstParameterProjectionPoint () const
//{
//  return _parametersProjectionPoint[0];
//}
//
//double FindClosestEdge:: getSecondParameterProjectionPoint () const
//{
//  return _parametersProjectionPoint[1];
//}

void FindClosestEdge::find(mesh::Edge &edge)
{
  PRECICE_TRACE(edge.vertex(0).getCoords(), edge.vertex(1).getCoords());

  // Methodology of book "Computational Geometry", Joseph O' Rourke, Chapter 7.2
  auto &a    = edge.vertex(0).getCoords();
  auto &b    = edge.vertex(1).getCoords();
  auto &norm = edge.getNormal();

  auto ret = math::barycenter::calcBarycentricCoordsForEdge(a, b, norm, _searchPoint);
  PRECICE_ASSERT(ret.barycentricCoords.size() == 2);

  bool inside = not(ret.barycentricCoords.array() < -math::NUMERICAL_ZERO_DIFFERENCE).any();

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
      _closestEdge                  = &edge;
    }
  }
}

} // namespace query
} // namespace precice
