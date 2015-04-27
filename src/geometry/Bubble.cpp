// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Bubble.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include <list>
#include "tarch/la/Scalar.h"

namespace precice {
namespace geometry {

tarch::logging::Log Bubble:: _log ( "precice::geometry::Bubble" );

Bubble:: Bubble
(
  const utils::DynVector& offset,
  double                  discretizationWidth,
  double                  radius,
  double                  deformation )
:
  Geometry ( offset ),
  _discretizationWidth ( discretizationWidth ),
  _radius ( radius ),
  _deformation ( deformation )
{
}

void Bubble:: specializedCreate
(
  mesh::Mesh& seed )
{
  preciceTrace ( "specializedCreate()" );
  using namespace mesh;
  int dimensions = seed.getDimensions();
  assertion1 ( (dimensions == 2) || (dimensions == 3), dimensions );
  if ( dimensions == 2 ){
    double currentRadius = _radius * (1.0 - _deformation / 4.0 + _deformation * 0.5
                           * (3.0 * std::pow(std::cos (0.0),2) - 1.0) );
    utils::Vector2D start(0.0);
    start(0) += currentRadius;
    Vertex* oldUpper = &seed.createVertex ( -1.0 * start );
    Vertex * oldLower = oldUpper;

    int numLatitudinal = getNumberLatitudinalElements (_discretizationWidth)+1;
    double angleInc = tarch::la::PI / (double) (numLatitudinal-1.0);

    double x, y;
    for (int i=1; i < numLatitudinal-1; i++) {
      currentRadius = _radius * (1.0 - _deformation / 4.0 + _deformation * 0.5
                           * (3.0 * std::pow(std::cos (i * angleInc),2) - 1.0) );
      x = - currentRadius * std::cos (i * angleInc);
      y =   currentRadius * std::sin (i * angleInc);
      Vertex* newUpper = & seed.createVertex ( utils::Vector2D(x, y) );
      Vertex* newLower = & seed.createVertex ( utils::Vector2D(x, -y) );
      seed.createEdge ( *newUpper, *oldUpper );
      seed.createEdge ( *oldLower, *newLower );
      oldUpper = newUpper;
      oldLower = newLower;
    }
    Vertex& last = seed.createVertex (start);
    seed.createEdge (last, *oldUpper);
    seed.createEdge (*oldLower, last);
  }
  else {
    preciceError ( "specializedCreate()", "Geometry bubble only available in 2D!");
  }
}

mesh::Vertex* Bubble:: getVertex
(
  mesh::Vertex&                               v0,
  mesh::Vertex&                               v1,
  std::map<std::pair<int,int>,mesh::Vertex*>& dividedEdges,
  mesh::Mesh&                                 seed )
{
  std::map<std::pair<int,int>,mesh::Vertex*>::iterator iter;
  std::pair<int,int> index = std::make_pair ( v0.getID(), v1.getID() );
  iter = dividedEdges.find ( index );
  if ( iter != dividedEdges.end() ) {
    return iter->second;
  }
  index = std::make_pair ( v1.getID(), v0.getID() );
  iter = dividedEdges.find ( index );
  if ( iter != dividedEdges.end() ) {
    return iter->second;
  }
  utils::Vector2D coords ( v0.getCoords() );
  coords += v1.getCoords();
  coords /= 2.0;
  coords *= _radius / tarch::la::norm2(coords);
  mesh::Vertex* vertex = & seed.createVertex ( coords );
  dividedEdges[index] = vertex;
  return vertex;
}

mesh::Edge* Bubble:: getEdge
(
  mesh::Vertex&                             v0,
  mesh::Vertex&                             v1,
  std::map<std::pair<int,int>,mesh::Edge*>& edges,
  mesh::Mesh&                               seed )
{
  std::map<std::pair<int,int>, mesh::Edge*>::iterator iter;
  std::pair<int,int> index = std::make_pair ( v0.getID(),v1.getID() );
  iter = edges.find ( index );
  if ( iter != edges.end() ) {
    return iter->second;
  }
  index = std::make_pair ( v1.getID(),v0.getID() );
  iter = edges.find ( index );
  if ( iter != edges.end() ) {
    return iter->second;
  }
  mesh::Edge* edge = & seed.createEdge ( v0, v1 );
  edges[index] = edge;
  return edge;
}

int Bubble:: getNumberLongitudinalElements
(
   double discretizationWidth ) const
{
  return 2; // Only for 2D
}

int Bubble:: getNumberLatitudinalElements
(
  const double discretizationWidth ) const
{
  return static_cast<int>((1.6 * _radius * tarch::la::PI / discretizationWidth) + 0.5);
}

}} // namespace precice, geometry
