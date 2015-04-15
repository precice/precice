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
#include <vector>

namespace precice {
namespace com {

tarch::logging::Log CommunicateMesh:: _log ( "precice::com::CommunicateMesh" );

CommunicateMesh:: CommunicateMesh
(
  com::Communication::SharedPointer communication )
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

  int numberOfVertices = mesh.vertices().size();
  _communication->send ( numberOfVertices, rankReceiver );
  if(numberOfVertices>0){
    double itemsToSend[numberOfVertices*dim];
    for(int i=0;i<numberOfVertices;i++){
      for(int d=0;d<dim;d++){
        itemsToSend[i*dim+d] = mesh.vertices()[i].getCoords()[d];
      }
    }
    _communication->send(itemsToSend,numberOfVertices*dim,rankReceiver);
  }

  int numberOfEdges = mesh.edges().size();
  _communication->send ( numberOfEdges, rankReceiver );
  if(numberOfEdges>0){
    //we need to send the vertexIDs first such that the right edges can be created later
    //contrary to the normal sendMesh, this variant must also work for adding delta meshes
    int vertexIDs[numberOfVertices];
    for(int i=0;i<numberOfVertices;i++){
      vertexIDs[i] = mesh.vertices()[i].getID();
    }
    _communication->send(vertexIDs,numberOfVertices,rankReceiver);

    int edgeIDs[numberOfEdges*2];
    for(int i=0;i<numberOfEdges;i++){
      edgeIDs[i*2] = mesh.edges()[i].vertex(0).getID();
      edgeIDs[i*2+1] = mesh.edges()[i].vertex(1).getID();
    }
    _communication->send(edgeIDs,numberOfEdges*2,rankReceiver);
  }

  if ( dim == 3 ) {
    int numberOfTriangles = mesh.triangles().size();
    _communication->send ( numberOfTriangles, rankReceiver );
    if(numberOfTriangles>0){
      //we need to send the edgeIDs first such that the right edges can be created later
      //contrary to the normal sendMesh, this variant must also work for adding delta meshes
      int edgeIDs[numberOfEdges];
      for(int i=0;i<numberOfEdges;i++){
        edgeIDs[i] = mesh.edges()[i].getID();
      }
      _communication->send(edgeIDs,numberOfEdges,rankReceiver);

      int triangleIDs[numberOfTriangles*3];
      for(int i=0;i<numberOfTriangles;i++){
        triangleIDs[i*3] = mesh.triangles()[i].edge(0).getID();
        triangleIDs[i*3+1] = mesh.triangles()[i].edge(1).getID();
        triangleIDs[i*3+2] = mesh.triangles()[i].edge(2).getID();
      }
      _communication->send(triangleIDs,numberOfTriangles*3,rankReceiver);
    }
  }
}

void CommunicateMesh:: receiveMesh
(
  mesh::Mesh& mesh,
  int         rankSender )
{
  preciceTrace2 ( "receiveMesh()", mesh.getName(), rankSender );
  using tarch::la::raw;
  int dim = mesh.getDimensions();

  std::vector<mesh::Vertex*> vertices;
  std::map<int, mesh::Vertex*> vertexMap;
  int numberOfVertices = 0;
  _communication->receive ( numberOfVertices, rankSender);

  if(numberOfVertices>0){
    double vertexCoords[numberOfVertices*dim];
    _communication->receive(vertexCoords,numberOfVertices*dim,rankSender);
    for ( int i=0; i < numberOfVertices; i++ ){
      utils::DynVector coords(dim);
      for ( int d=0; d < dim; d++){
        coords[d] = vertexCoords[i*dim+d];
      }
      mesh::Vertex& v = mesh.createVertex ( coords );
      assertion1 ( v.getID() >= 0, v.getID() );
      vertices.push_back(&v);
    }
  }

  int numberOfEdges = 0;
  std::vector<mesh::Edge*> edges;
  _communication->receive ( numberOfEdges, rankSender);
  if(numberOfEdges>0){
    int vertexIDs[numberOfVertices];
    _communication->receive(vertexIDs,numberOfVertices,rankSender);
    for( int i=0; i < numberOfVertices; i++){
      vertexMap[vertexIDs[i]] = vertices[i];
    }

    int edgeIDs[numberOfEdges];
    _communication->receive(edgeIDs,numberOfEdges*2,rankSender);
    for( int i=0; i < numberOfEdges; i++){
      assertion ( vertexMap.find(edgeIDs[i*2]) != vertexMap.end() );
      assertion ( vertexMap.find(edgeIDs[i*2+1]) != vertexMap.end() );
      assertion ( edgeIDs[i*2] != edgeIDs[i*2+1] );
      mesh::Edge& e = mesh.createEdge ( *vertexMap[edgeIDs[i*2]], *vertexMap[edgeIDs[i*2+1]] );
      edges.push_back(&e);
    }
  }

  if ( dim == 3 ){
    int numberOfTriangles = 0;
    _communication->receive ( numberOfTriangles, rankSender );
    if (numberOfTriangles > 0){
      assertion ( (edges.size() > 0) || (numberOfTriangles == 0) );
      int edgeIDs[numberOfEdges];
      _communication->receive(edgeIDs,numberOfEdges,rankSender);
      std::map<int,mesh::Edge*> edgeMap;
      for( int i=0; i < numberOfEdges; i++){
        edgeMap[edgeIDs[i]] = edges[i];
      }

      int triangleIDs[numberOfTriangles];
      _communication->receive(triangleIDs,numberOfTriangles*3,rankSender);

      for( int i=0; i < numberOfTriangles; i++){
        assertion ( edgeMap.find(triangleIDs[i*3]) != edgeMap.end() );
        assertion ( edgeMap.find(triangleIDs[i*3+1]) != edgeMap.end() );
        assertion ( edgeMap.find(triangleIDs[i*3+2]) != edgeMap.end() );
        assertion ( triangleIDs[i*3] != triangleIDs[i*3+1] );
        assertion ( triangleIDs[i*3+1] != triangleIDs[i*3+2] );
        assertion ( triangleIDs[i*3+2] != triangleIDs[i*3] );
        mesh.createTriangle ( *edgeMap[triangleIDs[i*3]], *edgeMap[triangleIDs[i*3+1]], *edgeMap[triangleIDs[i*3+2]] );
      }
    }
  }
}

void CommunicateMesh:: sendBoundingBox (
  const mesh::Mesh::BoundingBox & bb,
  int                rankReceiver ){
  preciceTrace1 ( "sendBoundingBox()", rankReceiver );
  int dim = bb.size();
  for(int d=0; d<dim; d++){
    _communication->send(bb[d].first, rankReceiver);
    _communication->send(bb[d].second, rankReceiver);
  }
}

void CommunicateMesh:: receiveBoundingBox (
  mesh::Mesh::BoundingBox & bb,
  int          rankSender ){
  preciceTrace1 ( "receiveBoundingBox()", rankSender );
  int dim = bb.size();
  for(int d=0; d<dim; d++){
    _communication->receive(bb[d].first, rankSender);
    _communication->receive(bb[d].second, rankSender);
  }
}

}} // namespace precice, com
