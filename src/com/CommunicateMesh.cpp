// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CommunicateMesh.hpp"
#include "Communication.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include <map>

namespace precice {
namespace com {

tarch::logging::Log CommunicateMesh:: _log ( "precice::com::CommunicateMesh" );

CommunicateMesh:: CommunicateMesh
(
  com::PtrCommunication communication )
:
  _communication ( communication )
{}

void CommunicateMesh:: sendMesh
(
  const mesh::Mesh& mesh,
  int               rankReceiver )
{
  preciceTrace2 ( "sendGeometry()", mesh.getName(), rankReceiver );
  using tarch::la::raw;
  int dim = mesh.getDimensions();
  _communication->startSendPackage ( rankReceiver );
  _communication->send ( (int)mesh.vertices().size(), rankReceiver );
  utils::DynVector coord(dim);
  foreach ( const mesh::Vertex& vertex, mesh.vertices() ){
    coord = vertex.getCoords();
    _communication->send (raw(coord), dim, rankReceiver);
  }
  _communication->send ( (int)mesh.edges().size(), rankReceiver );
  foreach ( const mesh::Edge& edge, mesh.edges() ){
    _communication->send ( edge.vertex(0).getID(), rankReceiver);
    _communication->send ( edge.vertex(1).getID(), rankReceiver);
  }
  if ( dim == 3 ) {
    _communication->send ( (int)mesh.triangles().size(), rankReceiver );
    foreach ( const mesh::Triangle& triangle, mesh.triangles() ){
      _communication->send ( triangle.edge(0).getID(), rankReceiver );
      _communication->send ( triangle.edge(1).getID(), rankReceiver );
      _communication->send ( triangle.edge(2).getID(), rankReceiver );
    }
  }
  _communication->finishSendPackage ();
}

void CommunicateMesh:: receiveMesh
(
  mesh::Mesh& mesh,
  int         rankSender )
{
  preciceTrace2 ( "receiveMesh()", mesh.getName(), rankSender );
  using tarch::la::raw;
  int dim = mesh.getDimensions();
  std::map<int, mesh::Vertex*> vertexMap;
  int vertexCount;
  _communication->startReceivePackage ( rankSender );
  _communication->receive ( vertexCount, rankSender);
  utils::DynVector coords(dim);
  for ( int i=0; i < vertexCount; i++ ){
    _communication->receive ( raw(coords), dim, rankSender );
    mesh::Vertex& v = mesh.createVertex ( coords );
    assertion1 ( v.getID() >= 0, v.getID() );
    vertexMap[v.getID()] = &v;
  }
  std::map<int,mesh::Edge*> edgeMap; // only used in 3D
  int edgeCount;
  _communication->receive ( edgeCount, rankSender );
  for ( int i=0; i < edgeCount; i++ ){
    int vertexIndex1 = -1;
    int vertexIndex2 = -1;
    _communication->receive ( vertexIndex1, rankSender );
    _communication->receive ( vertexIndex2, rankSender );

    assertion ( vertexMap.find(vertexIndex1) != vertexMap.end() );
    assertion ( vertexMap.find(vertexIndex2) != vertexMap.end() );
    assertion1 ( vertexIndex1 != vertexIndex2, vertexIndex1 );

    if ( dim == 2 ){
      mesh.createEdge ( *vertexMap[vertexIndex1], *vertexMap[vertexIndex2] );
    }
    else {
      assertion1 ( dim == 3, dim );
      mesh::Edge & e = mesh.createEdge (
          *vertexMap[vertexIndex1], *vertexMap[vertexIndex2] );
      assertion1 ( e.getID() >= 0, e.getID() );
      edgeMap[e.getID()] = &e;
    }
  }
  if ( dim == 3 ){
    int triangleCount = 0;
    assertion ( (edgeMap.size() > 0) || (triangleCount == 0) );
    _communication->receive ( triangleCount, rankSender );
    for ( int i=0; i < triangleCount; i++ ) {
      int index1 = -1;
      int index2 = -1;
      int index3 = -1;
      _communication->receive ( index1, rankSender );
      _communication->receive ( index2, rankSender );
      _communication->receive ( index3, rankSender );

      assertion ( edgeMap.find(index1) != edgeMap.end() );
      assertion ( edgeMap.find(index2) != edgeMap.end() );
      assertion ( edgeMap.find(index3) != edgeMap.end() );
      assertion ( index1 != index2 );
      assertion ( index2 != index3 );
      assertion ( index3 != index1 );

      mesh.createTriangle ( *edgeMap[index1], *edgeMap[index2], *edgeMap[index3] );
    }
  }
  _communication->finishReceivePackage ();
}

}} // namespace precice, com
