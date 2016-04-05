// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Quad.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "boost/assign.hpp"

namespace precice {
namespace mesh {

Quad:: Quad
(
  Edge& edgeOne,
  Edge& edgeTwo,
  Edge& edgeThree,
  Edge& edgeFour,
  int   id )
:
  PropertyContainer(),
  _edges( boost::assign::list_of(&edgeOne)(&edgeTwo)(&edgeThree)(&edgeFour).to_array(_edges) ),
  _vertexMap(),
  _id( id ),
  _normal( edgeOne.getDimensions() ),
  _center( edgeOne.getDimensions() ),
  _enclosingRadius ( 0.0 )
{
  assertion(edgeOne.getDimensions() == edgeTwo.getDimensions(),
             edgeOne.getDimensions(), edgeTwo.getDimensions() );
  assertion(edgeTwo.getDimensions() == edgeThree.getDimensions(),
             edgeTwo.getDimensions(), edgeThree.getDimensions() );
  assertion(edgeThree.getDimensions() == edgeFour.getDimensions(),
             edgeThree.getDimensions(), edgeFour.getDimensions() );
  assertion(getDimensions() == 3, getDimensions());

  // Determine vertex map
  Vertex& v0 = edge(0).vertex(0);
  Vertex& v1 = edge(0).vertex(1);

  // Check for edges 0 and 1 which vertex establishes connection
  if (&edge(1).vertex(0) == &v0){
    _vertexMap[0] = 1;
    _vertexMap[1] = 0;
  }
  else if (&edge(1).vertex(1) == &v0){
    _vertexMap[0] = 1;
    _vertexMap[1] = 1;
  }
  else if (&edge(1).vertex(0) == &v1){
    _vertexMap[0] = 0;
    _vertexMap[1] = 0;
  }
  else {
    assertion(&edge(1).vertex(1) == &v1);
    _vertexMap[0] = 0;
    _vertexMap[1] = 1;
  }

  // Check for edges 1 and 2 which vertex establishes connection
  if (_vertexMap[1] == 0){
    if (&edge(2).vertex(0) == &edge(1).vertex(1)){
      _vertexMap[2] = 0;
    }
    else {
      assertion(&edge(2).vertex(1) == &edge(1).vertex(1));
      _vertexMap[2] = 1;
    }
  }
  else if (_vertexMap[1] == 1){
    if (&edge(2).vertex(0) == &edge(1).vertex(0)){
      _vertexMap[2] = 0;
    }
    else {
      assertion(&edge(2).vertex(1) == &edge(1).vertex(0));
      _vertexMap[2] = 1;
    }
  }

  // Check for edges 2 and 3 which vertex establishes connection
  if (_vertexMap[2] == 0){
    if (&edge(3).vertex(0) == &edge(2).vertex(1)){
      _vertexMap[3] = 0;
    }
    else {
      assertion(&edge(3).vertex(1) == &edge(2).vertex(1));
      _vertexMap[3] = 1;
    }
  }
  else if (_vertexMap[2] == 1){
    if (&edge(3).vertex(0) == &edge(2).vertex(0)){
      _vertexMap[3] = 0;
    }
    else {
      assertion(&edge(3).vertex(1) == &edge(2).vertex(0));
      _vertexMap[3] = 1;
    }
  }

  assertion(&vertex(0) != &vertex(1));
  assertion(&vertex(0) != &vertex(2));
  assertion(&vertex(0) != &vertex(3));
  assertion(&vertex(1) != &vertex(2));
  assertion(&vertex(1) != &vertex(3));
  assertion(&vertex(2) != &vertex(3));

  assertion((_vertexMap[0] == 0) || (_vertexMap[0] == 1), _vertexMap[0]);
  assertion((_vertexMap[1] == 0) || (_vertexMap[1] == 1), _vertexMap[1]);
  assertion((_vertexMap[2] == 0) || (_vertexMap[2] == 1), _vertexMap[2]);
  assertion((_vertexMap[2] == 0) || (_vertexMap[2] == 1), _vertexMap[3]);
}

int Quad:: getDimensions() const
{
  return _edges[0]->getDimensions();
}

void Quad:: setEnclosingRadius
(
  double radius )
{
  _enclosingRadius = radius;
}

const utils::DynVector& Quad:: getNormal() const
{
  return _normal;
}

const utils::DynVector& Quad:: getCenter() const
{
  return _center;
}

double Quad:: getEnclosingRadius() const
{
  return _enclosingRadius;
}

}} // namespace precice, mesh
