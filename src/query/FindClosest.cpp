#include "FindClosest.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
#include "logging/LogMacros.hpp"
#include "math/barycenter.hpp"
#include "math/differences.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "query/FindClosestEdge.hpp"
#include "query/FindClosestTriangle.hpp"
#include "query/FindClosestVertex.hpp"

namespace precice {
namespace query {

std::ostream &operator<<(std::ostream &out, const InterpolationElement &val)
{
  out << '(' << *val.element << ", w:" << val.weight << ')';
  return out;
}

InterpolationElements generateInterpolationElements(
    const mesh::Vertex & /*location*/,
    const mesh::Vertex &element)
{
  return {{element, 1.0}};
}

InterpolationElements generateInterpolationElements(
    const mesh::Vertex &location,
    const mesh::Edge &  element)
{
  auto &A = element.vertex(0);
  auto &B = element.vertex(1);

  const auto bcoords = math::barycenter::calcBarycentricCoordsForEdge(
                           A.getCoords(),
                           B.getCoords(),
                           element.getNormal(),
                           location.getCoords())
                           .barycentricCoords;

  InterpolationElements elems;
  elems.emplace_back(A, bcoords(0));
  elems.emplace_back(B, bcoords(1));
  return elems;
}

InterpolationElements generateInterpolationElements(
    const mesh::Vertex &  location,
    const mesh::Triangle &element)
{
  auto &A = element.vertex(0);
  auto &B = element.vertex(1);
  auto &C = element.vertex(2);

  const auto bcoords = math::barycenter::calcBarycentricCoordsForTriangle(
                           A.getCoords(),
                           B.getCoords(),
                           C.getCoords(),
                           element.getNormal(),
                           location.getCoords())
                           .barycentricCoords;

  InterpolationElements elems;
  elems.emplace_back(A, bcoords(0));
  elems.emplace_back(B, bcoords(1));
  elems.emplace_back(C, bcoords(2));
  return elems;
}

bool FindClosest::hasFound() const
{
  return not _closest.meshIDs.empty();
}

const ClosestElement &FindClosest::getClosest()
{
  return _closest;
}

double FindClosest::getEuclidianDistance()
{
  return std::abs(_closest.distance);
}

const Eigen::VectorXd &FindClosest::getSearchPoint() const
{
  return _searchpoint;
}

bool FindClosest::determineClosest()
{
  PRECICE_TRACE(_searchpoint);
  using math::greater;
  _closest          = ClosestElement(_searchpoint.size());
  _closest.distance = std::numeric_limits<double>::max();
  int closestType   = -1;
  if (greater(_closest.distance, _findClosestVertex.getEuclidianDistance())) {
    _closest.distance = _findClosestVertex.getEuclidianDistance();
    closestType       = 0;
  }
  if (greater(_closest.distance, _findClosestEdge.getEuclidianDistance())) {
    _closest.distance = _findClosestEdge.getEuclidianDistance();
    closestType       = 1;
  }
  if (greater(_closest.distance, _findClosestTriangle.getEuclidianDistance())) {
    _closest.distance = _findClosestTriangle.getEuclidianDistance();
    closestType       = 2;
  }
  // Assign all properties to _closest
  Eigen::VectorXd normal = Eigen::VectorXd::Zero(_searchpoint.size());
  if (closestType == 0) { // Vertex
    mesh::Vertex &vertex     = _findClosestVertex.getClosestVertex();
    _closest.vectorToElement = vertex.getCoords() - _searchpoint;
    normal                   = vertex.getNormal();
    InterpolationElement element;
    element.element = &vertex;
    element.weight  = 1.0;
    _closest.interpolationElements.push_back(element);
  } else if (closestType == 1) { // Edge
    mesh::Edge &edge         = _findClosestEdge.getClosestEdge();
    _closest.vectorToElement = _findClosestEdge.getVectorToProjectionPoint();
    normal                   = edge.getNormal();
    InterpolationElement element0, element1;
    element0.element = &edge.vertex(0);
    element1.element = &edge.vertex(1);
    element0.weight  = _findClosestEdge.getProjectionPointParameter(0);
    element1.weight  = _findClosestEdge.getProjectionPointParameter(1);
    _closest.interpolationElements.push_back(element0);
    _closest.interpolationElements.push_back(element1);
  } else if (closestType == 2) { // Triangle
    mesh::Triangle &triangle = _findClosestTriangle.getClosestTriangle();
    _closest.vectorToElement = _findClosestTriangle.getVectorToProjectionPoint();
    normal                   = triangle.getNormal();
    InterpolationElement element0, element1, element2;
    element0.element = &triangle.vertex(0);
    element1.element = &triangle.vertex(1);
    element2.element = &triangle.vertex(2);
    element0.weight  = _findClosestTriangle.getProjectionPointParameter(0);
    element1.weight  = _findClosestTriangle.getProjectionPointParameter(1);
    element2.weight  = _findClosestTriangle.getProjectionPointParameter(2);
    _closest.interpolationElements.push_back(element0);
    _closest.interpolationElements.push_back(element1);
    _closest.interpolationElements.push_back(element2);
  } else {
    return false;
  }
  // Flip sign of distance, depending on normal of closest element
  if (_closest.vectorToElement.dot(normal) > 0.0) {
    _closest.distance *= -1.0;
  }
  return true;
}

void FindClosest::reset()
{
  _findClosestVertex   = FindClosestVertex(_searchpoint);
  _findClosestEdge     = FindClosestEdge(_searchpoint);
  _findClosestTriangle = FindClosestTriangle(_searchpoint);
}

} // namespace query
} // namespace precice
