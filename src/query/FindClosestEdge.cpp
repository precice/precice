// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "FindClosestEdge.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "utils/GeometryComputations.hpp"

namespace precice {
namespace query {

tarch::logging::Log FindClosestEdge:: _log ( "precice::query::FindClosestEdge" );

const utils::DynVector& FindClosestEdge:: getSearchPoint() const
{
  return _searchPoint;
}

bool FindClosestEdge:: hasFound() const
{
  return _closestEdge != NULL;
}

double FindClosestEdge:: getEuclidianDistance()
{
  return _shortestDistance;
}

mesh::Edge& FindClosestEdge:: getClosestEdge()
{
  assertion ( _closestEdge != NULL );
  return *_closestEdge;
}

const utils::DynVector& FindClosestEdge:: getVectorToProjectionPoint() const
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
  preciceTrace2 ( "find()", edge.vertex(0).getCoords(), edge.vertex(1).getCoords() );
  // Methodology of book "Computational Geometry", Joseph O' Rourke, Chapter 7.2
  boost::array<double,2> barycentricCoords;
  int dimensions = edge.getDimensions();
  assertion1 ( (dimensions == 2) || (dimensions == 3), dimensions );
  utils::DynVector projected(dimensions, 0.0);
  bool collinear = false;
  using utils::Vector2D; using utils::Vector3D;
  Vector2D a, b, ab, c, d;

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
    collinear = utils::GeometryComputations::collinear ( a, b, c );
    if ( collinear ) {
      // From p(s) = a + s(b-a) we get: s = (p(s) - a) / (b-a)
      int iMax = tarch::la::indexMax ( tarch::la::abs(ab) );
      assertion ( ! tarch::la::equals(ab(iMax), 0.0) );
      barycentricCoords[0] = (_searchPoint(iMax) - a(iMax)) / ab(iMax);
      barycentricCoords[1] = 1.0 - barycentricCoords[0];
      projected = _searchPoint;
    }
  }
  else { // 3D
    assertion1 ( dimensions == 3, dimensions );
    // Get parameters for parametric triangle representation: p(s) = a + s(b-a)
    Vector3D a3D = edge.vertex(0).getCoords();
    Vector3D b3D = edge.vertex(1).getCoords();
    Vector3D c3D = _searchPoint;
    Vector3D ab3D = b3D - a3D;
    Vector3D ac3D = c3D - a3D;
    collinear = utils::GeometryComputations::collinear ( a3D, b3D, c3D );
    if ( collinear ) {
      // From p(s) = a + s(b-a) we get: s = (p(s) - a) / (b-a)
      Vector3D ab3Dabs ( ab3D );
      int iMax = tarch::la::indexMax ( tarch::la::abs(ab3D, ab3Dabs) );
      assertion ( ! tarch::la::equals(ab3D[iMax], 0.0) );
      barycentricCoords[0] =
          (_searchPoint(iMax) - edge.vertex(0).getCoords()(iMax)) / ab3D(iMax);
      barycentricCoords[1] = 1.0 - barycentricCoords[0];
      projected = _searchPoint;
    }
    else {
      // Project parameters to 2D, where the projection plane is determined from
      // the normal direction, in order to prevent "faulty" projections.
      Vector3D normal;
      tarch::la::cross(ab3D,ac3D,normal);
      int indexToRemove = tarch::la::indexMax ( tarch::la::abs(normal,normal) );
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
        assertion1 ( indexToRemove == 2, indexToRemove );
        indices[0] = 0;
        indices[1] = 1;
      }
      assignList(a) = a3D[indices[0]], a3D[indices[1]];
      assignList(b) = b3D[indices[0]], b3D[indices[1]];
      assignList(ab) = ab3D[indices[0]], ab3D[indices[1]];
      assignList(c) = c3D[indices[0]], c3D[indices[1]];
      // 3D normal might be projected out, hence, compute new 2D edge normal
      Vector2D normal2D ( -1.0 * ab(1), ab(0) );
      d = c + normal2D;
    }
  }

  if ( not collinear ){
    // Compute denominator for solving 2x2 equation system
    double D = a(0)*(d(1)-c(1)) + b(0)*(c(1)-d(1)) + d(0)*ab(1) - c(0)*ab(1);
    assertion5 ( not tarch::la::equals(D, 0.0), a, b, c, d, ab );   // D == 0 would imply "normal // edge"

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
    if( barycentricCoords[i] < - tarch::la::NUMERICAL_ZERO_DIFFERENCE) {
      inside = false;
    }
  }

  // if valid, compute distance to triangle and evtl. store distance
  if (inside) {
    utils::DynVector distanceVector = projected;
    distanceVector -= _searchPoint;
    double distance = tarch::la::norm2 ( distanceVector );
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
