// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Merge.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace mesh {

Merge:: Merge()
:
  //PROPERTY_ID ( PropertyContainer::getFreePropertyID() ),
  _merged()
{}

Merge:: ~Merge()
{
//  foreach ( Vertex& vertex, _merged.vertices() ){
//    vertex.deleteProperty ( PROPERTY_ID );
//  }
//  foreach ( Edge& edge, _merged.edges() ){
//    edge.deleteProperty ( PROPERTY_ID );
//  }
//  foreach ( Triangle& triangle, _merged.triangles() ){
//    triangle.deleteProperty ( PROPERTY_ID );
//  }
}

Group& Merge:: content()
{
  return _merged;
}

//int Merge:: getPropertyIndex()
//{
//  return PROPERTY_ID;
//}
//
//bool Merge:: isUnique
//(
//  PropertyContainer& cont )
//{
//  if ( cont.hasProperty(PROPERTY_ID) ){
//    return false;
//  }
//  cont.setProperty ( PROPERTY_ID, 1 );
//  return true;
//}

}} // namespace precice, mesh
