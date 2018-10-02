#include "FindClosest.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Quad.hpp"
#include "math/math.hpp"
#include <limits>

namespace precice {
namespace query {

bool FindClosest:: hasFound() const
{
  return not _closest.meshIDs.empty();
}

const ClosestElement& FindClosest:: getClosest()
{
  return _closest;
}

double FindClosest:: getEuclidianDistance()
{
  return std::abs(_closest.distance);
}

const Eigen::VectorXd& FindClosest:: getSearchPoint() const
{
  return _searchpoint;
}

bool FindClosest:: determineClosest()
{
  TRACE(_searchpoint);
  using math::greater;
  _closest = ClosestElement(_searchpoint.size());
  _closest.distance = std::numeric_limits<double>::max();
  int closestType = -1;
  if ( greater(_closest.distance, _findClosestVertex.getEuclidianDistance()) ) {
    _closest.distance = _findClosestVertex.getEuclidianDistance();
    closestType = 0;
  }
  if ( greater(_closest.distance, _findClosestEdge.getEuclidianDistance()) ) {
    _closest.distance = _findClosestEdge.getEuclidianDistance();
    closestType = 1;
  }
  if ( greater(_closest.distance, _findClosestTriangle.getEuclidianDistance()) ) {
    _closest.distance = _findClosestTriangle.getEuclidianDistance();
    closestType = 2;
  }
  if ( greater(_closest.distance, _findClosestQuad.getEuclidianDistance()) ) {
    _closest.distance = _findClosestQuad.getEuclidianDistance();
    closestType = 3;
  }
  // Assign all properties to _closest
  Eigen::VectorXd normal = Eigen::VectorXd::Zero(_searchpoint.size());
  if ( closestType == 0 ) { // Vertex
    mesh::Vertex& vertex = _findClosestVertex.getClosestVertex ();
    vertex.getProperties ( vertex.INDEX_GEOMETRY_ID, _closest.meshIDs );
    _closest.vectorToElement = vertex.getCoords() - _searchpoint;
    normal = vertex.getNormal();
    InterpolationElement element;
    element.element = & vertex;
    element.weight = 1.0;
    _closest.interpolationElements.push_back ( element );
  }
  else if ( closestType == 1 ) { // Edge
    mesh::Edge& edge = _findClosestEdge.getClosestEdge();
    edge.getProperties ( edge.INDEX_GEOMETRY_ID, _closest.meshIDs );
    _closest.vectorToElement = _findClosestEdge.getVectorToProjectionPoint();
    normal = edge.getNormal();
    InterpolationElement element0, element1;
    element0.element = & edge.vertex (0);
    element1.element = & edge.vertex (1);
    element0.weight = _findClosestEdge.getProjectionPointParameter(0);
    element1.weight = _findClosestEdge.getProjectionPointParameter(1);
    _closest.interpolationElements.push_back ( element0 );
    _closest.interpolationElements.push_back ( element1 );
  }
  else if ( closestType == 2 ) { // Triangle
    mesh::Triangle& triangle = _findClosestTriangle.getClosestTriangle ();
    triangle.getProperties ( triangle.INDEX_GEOMETRY_ID, _closest.meshIDs );
    _closest.vectorToElement = _findClosestTriangle.getVectorToProjectionPoint();
    normal = triangle.getNormal();
    InterpolationElement element0, element1, element2;
    element0.element = & triangle.vertex (0);
    element1.element = & triangle.vertex (1);
    element2.element = & triangle.vertex (2);
    element0.weight = _findClosestTriangle.getProjectionPointParameter(0);
    element1.weight = _findClosestTriangle.getProjectionPointParameter(1);
    element2.weight = _findClosestTriangle.getProjectionPointParameter(2);
    _closest.interpolationElements.push_back ( element0 );
    _closest.interpolationElements.push_back ( element1 );
    _closest.interpolationElements.push_back ( element2 );
  }
  else if ( closestType == 3 ) { // Quad
    mesh::Quad& quad = _findClosestQuad.getClosestQuad();
    quad.getProperties(quad.INDEX_GEOMETRY_ID, _closest.meshIDs);
    _closest.vectorToElement = _findClosestQuad.getVectorToProjectionPoint();
    normal = quad.getNormal();
    InterpolationElement element0, element1, element2, element3;
    element0.element = &quad.vertex(0);
    element1.element = &quad.vertex(1);
    element2.element = &quad.vertex(2);
    element3.element = &quad.vertex(3);
    element0.weight = _findClosestQuad.getProjectionPointParameter(0);
    element1.weight = _findClosestQuad.getProjectionPointParameter(1);
    element2.weight = _findClosestQuad.getProjectionPointParameter(2);
    element3.weight = _findClosestQuad.getProjectionPointParameter(3);
    _closest.interpolationElements.push_back(element0);
    _closest.interpolationElements.push_back(element1);
    _closest.interpolationElements.push_back(element2);
    _closest.interpolationElements.push_back(element3);
  }
  else {
    return false;
  }
// Flip sign of distance, depending on normal of closest element
if (_closest.vectorToElement.dot(normal) > 0.0) {
    _closest.distance *= -1.0;
  }
  return true;
}

void FindClosest:: reset()
{
  _findClosestVertex = FindClosestVertex(_searchpoint);
  _findClosestEdge = FindClosestEdge(_searchpoint);
  _findClosestTriangle = FindClosestTriangle(_searchpoint);
  _findClosestQuad = FindClosestQuad(_searchpoint);
}

}} // namespace precice, query
