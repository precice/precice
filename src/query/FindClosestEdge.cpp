#include "FindClosestEdge.hpp"
#include <Eigen/Dense>
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "math/geometry.hpp"

namespace precice {
namespace query {

FindClosestEdge:: FindClosestEdge
(
  const Eigen::VectorXd& searchPoint )
:
  _searchPoint ( searchPoint ),
  _vectorToProjectionPoint ( Eigen::VectorXd::Constant(searchPoint.size(), std::numeric_limits<double>::max()) ),
  _parametersProjectionPoint( {_shortestDistance , _shortestDistance} )
{
  assertion ( (_searchPoint.size() == 2) || (_searchPoint.size() == 3),
               _searchPoint.size() );
}

const Eigen::VectorXd& FindClosestEdge:: getSearchPoint() const
{
  return _searchPoint;
}

bool FindClosestEdge:: hasFound() const
{
  return _closestEdge != nullptr;
}

double FindClosestEdge:: getEuclidianDistance()
{
  return _shortestDistance;
}

mesh::Edge& FindClosestEdge:: getClosestEdge()
{
  assertion ( _closestEdge != nullptr );
  return *_closestEdge;
}

const Eigen::VectorXd& FindClosestEdge:: getVectorToProjectionPoint() const
{
  return _vectorToProjectionPoint;
}

double FindClosestEdge:: getProjectionPointParameter
(
  int index ) const
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

void FindClosestEdge:: find ( mesh::Edge& edge )
{
  TRACE(edge.vertex(0).getCoords(), edge.vertex(1).getCoords() );
  using Eigen::Vector2d; using Eigen::Vector3d;
  // Methodology of book "Computational Geometry", Joseph O' Rourke, Chapter 7.2
  std::array<double,2> barycentricCoords;
  int dimensions = edge.getDimensions();
  assertion ( (dimensions == 2) || (dimensions == 3), dimensions );
  Eigen::VectorXd projected = Eigen::VectorXd::Zero(dimensions);
  bool collinear = false;
  Vector2d a, b, ab, c, d;

  if ( dimensions == 2 ){
    // Get parameters for parametric edge representation: p(s) = a + s(b-a)
    a = edge.vertex(0).getCoords();
    b = edge.vertex(1).getCoords();
    ab = b;
    ab -= a;
    // Same for intersecting normal from searchpoint: q(t) = c + t(d - c)
    c = _searchPoint;
    d = _searchPoint;
    d += edge.getNormal();
    collinear = math::geometry::collinear ( a, b, c );
    if ( collinear ) {
      // From p(s) = a + s(b-a) we get: s = (p(s) - a) / (b-a)
      int iMax;
      ab.cwiseAbs().maxCoeff(&iMax);
      assertion ( ! math::equals(ab(iMax), 0.0) );
      barycentricCoords[0] = (_searchPoint(iMax) - a(iMax)) / ab(iMax);
      barycentricCoords[1] = 1.0 - barycentricCoords[0];
      projected = _searchPoint;
    }
  }
  else { // 3D
    assertion ( dimensions == 3, dimensions );
    // Get parameters for parametric triangle representation: p(s) = a + s(b-a)
    Vector3d a3D = edge.vertex(0).getCoords();
    Vector3d b3D = edge.vertex(1).getCoords();
    Vector3d c3D = _searchPoint;
    Vector3d ab3D = b3D - a3D;
    Vector3d ac3D = c3D - a3D;
    collinear = math::geometry::collinear ( a3D, b3D, c3D );
    if ( collinear ) {
      // From p(s) = a + s(b-a) we get: s = (p(s) - a) / (b-a)
      int iMax;
      ab3D.cwiseAbs().maxCoeff(&iMax);
      assertion ( ! math::equals(ab3D[iMax], 0.0) );
      barycentricCoords[0] =
        (_searchPoint(iMax) - edge.vertex(0).getCoords()(iMax)) / ab3D(iMax);
      barycentricCoords[1] = 1.0 - barycentricCoords[0];
      projected = _searchPoint;
    }
    else {
      // Project parameters to 2D, where the projection plane is determined from
      // the normal direction, in order to prevent "faulty" projections.
      Vector3d normal;
      normal = ab3D.cross(ac3D);
      int indexToRemove;
      normal.cwiseAbs().maxCoeff(&indexToRemove);
      int indices[2];
      if (indexToRemove == 0){
        indices[0] = 1;
        indices[1] = 2;
      }
      else if (indexToRemove == 1){
        indices[0] = 0;
        indices[1] = 2;
      }
      else {
        assertion ( indexToRemove == 2, indexToRemove );
        indices[0] = 0;
        indices[1] = 1;
      }
      a << a3D[indices[0]], a3D[indices[1]];
      b << b3D[indices[0]], b3D[indices[1]];
      ab << ab3D[indices[0]], ab3D[indices[1]];
      c << c3D[indices[0]], c3D[indices[1]];
      // 3D normal might be projected out, hence, compute new 2D edge normal
      Vector2d normal2D ( -1.0 * ab(1), ab(0) );
      d = c + normal2D;
    }
  }

  if ( not collinear ){
    // Compute denominator for solving 2x2 equation system
    double D = a(0)*(d(1)-c(1)) + b(0)*(c(1)-d(1)) + d(0)*ab(1) - c(0)*ab(1);
    assertion ( not math::equals(D, 0.0), a, b, c, d, ab );   // D == 0 would imply "normal // edge"

    // Compute triangle segment parameter s, which is in [0, 1], if the
    // intersection point is within the triangle.
    barycentricCoords[0] = (a(0)*(d(1)-c(1)) +
                           c(0)*(a(1)-d(1)) +
                           d(0)*(c(1)-a(1))) / D;
    barycentricCoords[1] = 1.0 - barycentricCoords[0];

    // Compute coordinates of projected point, which has to be done dimension
    // dependent again:
    if ( dimensions == 2 ){
      projected = ab;                    // = b - a
      projected *= barycentricCoords[0]; // = bary0 * (b - a)
      projected += a;                    // = a + bary0 * (b - a)
    }
    else {
      projected = edge.vertex(1).getCoords();  // = b
      projected -= edge.vertex(0).getCoords(); // = b - a
      projected *= barycentricCoords[0];       // = bary0 * (b - a)
      projected += edge.vertex(0).getCoords(); // = a + bary0 * (b - a)
    }
  }

  bool inside = true;
  for (int i=0; i < 2; i++) {
    if( barycentricCoords[i] < - math::NUMERICAL_ZERO_DIFFERENCE) {
      inside = false;
    }
  }

  // if valid, compute distance to triangle and evtl. store distance
  if (inside) {
    Eigen::VectorXd distanceVector = projected;
    distanceVector -= _searchPoint;
    double distance = distanceVector.norm();
    if ( _shortestDistance > distance ) {
      _shortestDistance = distance;
      _vectorToProjectionPoint = distanceVector;
      _parametersProjectionPoint[0] = barycentricCoords[1];
      _parametersProjectionPoint[1] = barycentricCoords[0];
      _closestEdge = &edge;
    }
  }
}

}} // namespace precice, query
