#include "FindClosestTriangle.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
#include "math/GeometryComputations.hpp"
#include "math/math.hpp"

namespace precice {
namespace query {

logging::Logger FindClosestTriangle:: _log ( "precice::query::FindClosestTriangle" );

FindClosestTriangle:: FindClosestTriangle
(
  const Eigen::VectorXd& searchPoint )
:
  _searchPoint ( searchPoint ),
  _shortestDistance ( std::numeric_limits<double>::max() ),
  _vectorToProjectionPoint ( Eigen::VectorXd::Constant(_searchPoint.size(), std::numeric_limits<double>::max()) ),
  _parametersProjectionPoint( {_shortestDistance, _shortestDistance, _shortestDistance } ),
  _closestTriangle ( nullptr )
{}

const Eigen::VectorXd& FindClosestTriangle:: getSearchPoint() const
{
  return _searchPoint;
}

bool FindClosestTriangle:: hasFound() const
{
  return _closestTriangle != nullptr;
}

double FindClosestTriangle:: getEuclidianDistance()
{
  return _shortestDistance;
}

mesh::Triangle & FindClosestTriangle:: getClosestTriangle()
{
  assertion(_closestTriangle != nullptr);
  return *_closestTriangle;
}

const Eigen::VectorXd& FindClosestTriangle:: getVectorToProjectionPoint() const
{
  return _vectorToProjectionPoint;
}

double FindClosestTriangle:: getProjectionPointParameter
(
  int index ) const
{
  return _parametersProjectionPoint[index];
}

void FindClosestTriangle:: find
(
  mesh::Triangle& triangle )
{
  using Eigen::Vector2d; using Eigen::Vector3d;
  
  Vector3d projected;

  // Methodology of book "Computational Geometry", Joseph O' Rourke, Chapter 7.3
  // with the barycentric coordinates method and real projection into 2D, instead
  // of outprojecting one coordinate

  // Parametric representation for triangle plane:
  // (x, y, z) * normal = d
  Vector3d a = triangle.vertex(0).getCoords();
  Vector3d normal = triangle.getNormal();
  double d = normal.dot(a);

  // Parametric description of line from searchpoint orthogonal to triangle:
  // _searchPoint + t * normal = x     (where t is parameter)
  // Determine t such that x lies on triangle plane:
  double t = d - _searchPoint.dot(normal) / normal.dot(normal);

  // Compute projected point with parameter t:
  projected = normal;
  projected *= t;
  projected += _searchPoint;

  // Project everything to 2D
  Vector3d normalAbs;
  int iMax;
  normal.cwiseAbs().maxCoeff(&iMax);
  int indices[2];
  if (iMax == 0){
    indices[0] = 1;
    indices[1] = 2;
  }
  else if (iMax == 1){
    indices[0] = 0;
    indices[1] = 2;
  }
  else {
    assertion(iMax == 2, iMax);
    indices[0] = 0;
    indices[1] = 1;
  }
  Vector2d a2D(triangle.vertex(0).getCoords()[indices[0]],
               triangle.vertex(0).getCoords()[indices[1]]);
  Vector2d b2D(triangle.vertex(1).getCoords()[indices[0]],
               triangle.vertex(1).getCoords()[indices[1]]);
  Vector2d c2D(triangle.vertex(2).getCoords()[indices[0]],
               triangle.vertex(2).getCoords()[indices[1]]);
  Vector2d projected2D(projected[indices[0]], projected[indices[1]]);
  // Compute barycentric coordinates by solving linear 3x3 system
  Vector3d rhs (projected2D(0), projected2D(1), 1);
  Eigen::Matrix<double, 3,3> A;
  A << a2D(0), b2D(0), c2D(0),
    a2D(1), b2D(1), c2D(1),
    1,      1,      1;
  Eigen::Vector3d barycentricCoords = A.colPivHouseholderQr().solve(rhs);
  // Determine from barycentric coordinates, if point is inside triangle
  bool inside = not (barycentricCoords.array() < - math::NUMERICAL_ZERO_DIFFERENCE).any();

  // If inside, compute distance to triangle and evtl. store distance
  if (inside) {
    Vector3d distanceVector = projected - _searchPoint;
    double distance = distanceVector.norm();
    if (_shortestDistance > distance){
      _shortestDistance = distance;
      _vectorToProjectionPoint = distanceVector;
      _parametersProjectionPoint[0] = barycentricCoords(0);
      _parametersProjectionPoint[1] = barycentricCoords(1);
      _parametersProjectionPoint[2] = barycentricCoords(2);
      _closestTriangle = &triangle;
    }
  }
}

}} // namespace precice, query
