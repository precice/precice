// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "FindClosestTriangle.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Dimensions.hpp"
#include "utils/GeometryComputations.hpp"
#include "tarch/la/Matrix.h"

namespace precice {
namespace query {

tarch::logging::Log FindClosestTriangle:: _log ( "precice::query::FindClosestTriangle" );

const utils::DynVector& FindClosestTriangle:: getSearchPoint() const
{
  return _searchPoint;
}

bool FindClosestTriangle:: hasFound() const
{
  return _closestTriangle != NULL;
}

double FindClosestTriangle:: getEuclidianDistance()
{
  return _shortestDistance;
}

mesh::Triangle & FindClosestTriangle:: getClosestTriangle()
{
  assertion(_closestTriangle != NULL);
  return *_closestTriangle;
}

const utils::DynVector& FindClosestTriangle:: getVectorToProjectionPoint() const
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
  typedef utils::GeometryComputations GeoComp;
  using utils::Vector2D; using utils::Vector3D;
  using namespace tarch::la;

  Vector3D barycentricCoords;
  Vector3D projected;

  // Methodology of book "Computational Geometry", Joseph O' Rourke, Chapter 7.3
  // with the barycentric coordinates method and real projection into 2D, instead
  // of outprojecting one coordinate

  // Parametric representation for triangle plane:
  // (x, y, z) * normal = d
  Vector3D a = triangle.vertex(0).getCoords();
  Vector3D normal = triangle.getNormal();
  double d = dot(normal, a);

  // Parametric description of line from searchpoint orthogonal to triangle:
  // _searchPoint + t * normal = x     (where t is parameter)
  // Determine t such that x lies on triangle plane:
  double t = d - dot(_searchPoint, normal) / dot(normal, normal);

  // Compute projected point with parameter t:
  projected = normal;
  projected *= t;
  projected += _searchPoint;

  // Project everything to 2D
  Vector3D normalAbs;
  int iMax = indexMax(abs(normal, normalAbs));
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
    assertion1(iMax == 2, iMax);
    indices[0] = 0;
    indices[1] = 1;
  }
  Vector2D a2D(triangle.vertex(0).getCoords()[indices[0]],
               triangle.vertex(0).getCoords()[indices[1]]);
  Vector2D b2D(triangle.vertex(1).getCoords()[indices[0]],
               triangle.vertex(1).getCoords()[indices[1]]);
  Vector2D c2D(triangle.vertex(2).getCoords()[indices[0]],
               triangle.vertex(2).getCoords()[indices[1]]);
  Vector2D projected2D(projected[indices[0]], projected[indices[1]]);
  // Compute barycentric coordinates by solving linear 3x3 system
  Vector3D rhs (projected2D(0), projected2D(1), 1);
  Matrix<3,3,double> A;
  assignList(A) = a2D(0), b2D(0), c2D(0),
    a2D(1), b2D(1), c2D(1),
    1,      1,      1;
  solveSystem3x3(A, rhs, barycentricCoords);

  // Determine from barycentric coordinates, if point is inside triangle
  bool inside = true;
  for (int i=0; i < 3; i++){
    if (barycentricCoords(i) < - NUMERICAL_ZERO_DIFFERENCE){
      inside = false;
    }
  }

  // if inside, compute distance to triangle and evtl. store distance
  if (inside){
    Vector3D distanceVector = projected - _searchPoint;
    double distance = norm2(distanceVector);
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
