// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "BalanceVertexPositionAction.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/Scalar.h"
#include "mesh/Mesh.hpp"
#include "mesh/Edge.hpp"
#include "utils/Dimensions.hpp"

namespace precice {
namespace impl {

tarch::logging::Log BalanceVertexPositionAction::
  _log ( "precice::impl::BalanceVertexPositionAction" );

BalanceVertexPositionAction:: BalanceVertexPositionAction
(
  TimingConstants       timing,
  const mesh::PtrMesh & mesh,
  double                convergenceLimit,
  int                   maxIterations )
:
  AbstractDataAction ( timing ),
  _mesh ( mesh ),
  _eps ( convergenceLimit ),
  _maxIterations ( maxIterations )
{
  assertion ( _mesh.use_count() > 0 );
  preciceCheck ( tarch::la::greater(_eps,0.0), "BalanceVertexPositionAction()",
      "Convergence limit has to be larger than " << tarch::la::NUMERICAL_ZERO_DIFFERENCE
      << " for balance-vertex-position action!");
  preciceCheck ( maxIterations > 0, "BalanceVertexPositionAction()",
      "Maximum iteration number has to be larger than 0 for "
      << "balance-vertex-position action!" );
}


void BalanceVertexPositionAction:: performAction
(
  double dt,
  double computedPartFullDt,
  double fullDt )
{
  preciceTrace ( "performAction()" );
  assertion ( _mesh.use_count() > 0 );
  tarch::la::DynamicVector<Vector> edgeVectors ( _mesh->vertices().size() );
  double errorMeasure = _eps * 10.0;
  int iterations = 0;
  while ( (errorMeasure > _eps) && (iterations < _maxIterations) ){
    errorMeasure = 0.0;
    tarch::la::assign(edgeVectors) = Vector(0.0);
    foreach ( mesh::Edge & edge, _mesh->edges() ){
      mesh::Vertex & v0 = edge.vertex(0);
      mesh::Vertex & v1 = edge.vertex(1);
      Vector ab ( v0.getCoords() );
      ab -= v1.getCoords();
      ab *= 0.5;
      assertion1 ( v0.getID() < edgeVectors.size(), v0.getID() );
      assertion1 ( v1.getID() < edgeVectors.size(), v1.getID() );
      //precicePrint ( "Addign " << ab << " to vertex " << v0.getID()  );
      //precicePrint ( "Subtracting " << ab << " from vertex " << v1.getID()  );
      edgeVectors[v0.getID()] -= ab;
      edgeVectors[v1.getID()] += ab;
    }
    foreach ( mesh::Vertex & vertex, _mesh->vertices() ){
      Vector point ( vertex.getCoords() );
      assertion1 ( vertex.getID() < edgeVectors.size(), vertex.getID() );
      point += edgeVectors[vertex.getID()];
      //precicePrint ( "Vertex " << vertex.getID() << " has edgeVectors = "
      //               << edgeVectors[vertex.getID()] );

#     ifdef Dim2
      // Get parameters for parametric representation of line orthogonal to vertex
      // normal: p(s) = a + s(b-a)
      Vector a = vertex.getCoords();
      Vector ab ( -vertex.getNormal()[1], vertex.getNormal()[0] );
      Vector b = a;
      b += ab;
      // Same for intersecting normal from point to line: q(t) = c + t(d - c)
      Vector c = point;
      Vector d = point;
      d += vertex.getNormal();
      //precicePrint ( "a = " << a << ", b = " << b << ", ab = " << ab );
      //precicePrint ( "c = " << c << ", d = " << d );
      // Compute denominator for solving 2x2 equation system
      double D = a(0)*(d(1)-c(1)) + b(0)*(c(1)-d(1)) + d(0)*ab(1) - c(0)*ab(1);
      assertion5 ( not tarch::la::equals(D, 0.0), a, b, c, d, ab );   // D == 0 would imply "normal // edge"
      // Compute parameter of intersection point on line ab
      double param = (a(0)*(d(1)-c(1)) + c(0)*(a(1)-d(1)) + d(0)*(c(1)-a(1))) / D;
      // Compute coordinates of projected point:
      Vector projectedPoint = ab;
      projectedPoint *= param;
      projectedPoint += a;
      //precicePrint ( "projectedPoint = " << projectedPoint );
#     else
      preciceError ( "performAction()",
          "BalanceVertexPositionAction is only supported for 2D!" );
#     endif
      Vector coordDelta ( projectedPoint );
      coordDelta -= vertex.getCoords();
      errorMeasure += tarch::la::dot(coordDelta, coordDelta);
      vertex.setCoords ( projectedPoint );
    }
    _mesh->computeState ();
    //precicePrint ( "Error measure = " << errorMeasure );
    errorMeasure = std::sqrt ( errorMeasure );
    iterations ++;
  }
}

}} // namespace precice, impl
