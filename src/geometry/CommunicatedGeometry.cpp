// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CommunicatedGeometry.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace geometry {

tarch::logging::Log CommunicatedGeometry:: _log ( "precice::geometry::CommunicatedGeometry" );

CommunicatedGeometry:: CommunicatedGeometry
(
  const utils::DynVector&  offset,
  const std::string&       accessor,
  const std::string&       provider,
  com::PtrCommunication   masterSlaveCom,
  int                     rank,
  int                     size,
  int                     dimensions)
:
  Geometry ( offset ),
  _accessorName ( accessor ),
  _providerName ( provider ),
  _receivers (),
  _masterSlaveCom(masterSlaveCom),
  _rank(rank),
  _size(size),
  _dimensions(dimensions),
  _boundingFromMapping(),
  _boundingToMapping()
{
  preciceTrace2 ( "CommunicatedGeometry()", accessor, provider );
}

void CommunicatedGeometry:: addReceiver
(
  const std::string&     receiver,
  com::PtrCommunication com )
{
  preciceTrace1 ( "addReceiver()", receiver );
  assertion ( com.get() != NULL );
  preciceCheck ( ! utils::contained(receiver, _receivers),
                 "addReceiver()", "Receiver \"" << receiver
                 << "\" has been added already to communicated geometry!" );
  preciceCheck ( receiver != _providerName, "addReceiver()",
                 "Receiver \"" << receiver << "\" cannot be the same as "
                 << "provider in communicated geometry!" );
  _receivers[receiver] = com;
}

void CommunicatedGeometry:: specializedCreate
(
  mesh::Mesh& seed )
{
  preciceTrace1 ( "specializedCreate()", seed.getName() );
  preciceCheck ( _receivers.size() > 0, "specializedCreate()",
                 "No receivers specified for communicated geometry to create "
                 << "mesh \"" << seed.getName() << "\"!" );
  if ( _accessorName == _providerName ) {
    if(_size>1){
      gatherMesh(seed);
    }
    if(_rank==0){
      preciceCheck ( seed.vertices().size() > 0,
                     "specializedCreate()", "Participant \"" << _accessorName
                     << "\" provides an invalid (possibly empty) mesh \""
                     << seed.getName() << "\"!" );
      typedef std::map<std::string,com::PtrCommunication>::value_type Pair;
      foreach ( Pair & pair, _receivers ) {
        if ( ! pair.second->isConnected() ) {
          pair.second->acceptConnection ( _providerName, pair.first, 0, 1 );
        }
        com::CommunicateMesh(pair.second).sendMesh ( seed, 0 );
      }
    }
  }
  else if ( utils::contained(_accessorName, _receivers) ) {
    if(_rank==0){
      assertion ( seed.vertices().size() == 0 );
      assertion ( utils::contained(_accessorName, _receivers) );
      com::PtrCommunication com ( _receivers[_accessorName] );
      if ( ! com->isConnected() ) {
        com->requestConnection ( _providerName, _accessorName, 0, 1 );
      }
      com::CommunicateMesh(com).receiveMesh ( seed, 0 );
    }
    if(_size>1){
      scatterMesh(seed);
    }
  }
  else {
    preciceError( "specializedCreate()", "Participant \"" << _accessorName
                  << "\" uses a communicated geometry to create mesh \""
                  << seed.getName()
                  << "\" but is neither provider nor receiver!" );
  }
}

void CommunicatedGeometry:: gatherMesh(
    mesh::Mesh& seed)
{
  preciceTrace1 ( "gatherMesh()", _rank );

  if(_rank>0){ //slave
    com::CommunicateMesh(_masterSlaveCom).sendMesh ( seed, 0 );
  }
  else{ //master
    assertion(_rank==0);
    assertion(_size>1);

    int numberOfVertices = 0;
    //vertices of master mesh part do already exist
    foreach ( const mesh::Vertex& vertex, seed.vertices() ){
      seed.getVertexDistribution()[0].push_back(numberOfVertices);
      numberOfVertices++;
    }

    mesh::Mesh slaveMesh("SlaveMesh", _dimensions, seed.isFlipNormals());
    mesh::Mesh& rSlaveMesh = slaveMesh;

    for(int rankSlave = 1; rankSlave < _size; rankSlave++){
      //slaves have ranks from 0 to size-2
      //TODO better rewrite accept/request connection
      rSlaveMesh.clear();
      com::CommunicateMesh(_masterSlaveCom).receiveMesh ( rSlaveMesh, rankSlave-1);
      utils::DynVector coord(_dimensions);
      seed.addMesh(rSlaveMesh);

      for(int i = 0; i < rSlaveMesh.vertices().size(); i++){
        seed.getVertexDistribution()[rankSlave].push_back(numberOfVertices);
        numberOfVertices++;
      }
    }
  }
}

void CommunicatedGeometry:: scatterMesh(
    mesh::Mesh& seed)
{
  preciceTrace1 ( "scatterMesh()", _rank );
  using tarch::la::raw;

  if(_rank>0){ //slave

    com::CommunicateMesh(_masterSlaveCom).receiveMesh ( seed, 0);

    seed.computeState();
    computeBoundingMappings();

    mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());

    preciceDebug("Bounding mesh. #vertices: " << seed.vertices().size()
         <<", #edges: " << seed.edges().size()
         <<", #triangles: " << seed.triangles().size() << ", rank: " << _rank);

    std::vector<int> globalVertexIDs;
    std::map<int, mesh::Vertex*> vertexMap;
    std::map<int, mesh::Edge*> edgeMap;

    foreach ( const mesh::Vertex& vertex, seed.vertices() ){
      if(doesVertexContribute(vertex.getID())){
        mesh::Vertex& v = filteredMesh.createVertex(vertex.getCoords());
        globalVertexIDs.push_back(vertex.getID());
        vertexMap[vertex.getID()] = &v;
      }
    }

    //add all edges formed by contributing vertices
    foreach ( mesh::Edge& edge, seed.edges() ){
      int vertexIndex1 = edge.vertex(0).getID();
      int vertexIndex2 = edge.vertex(1).getID();
      if(vertexMap.find(vertexIndex1) != vertexMap.end() &&
         vertexMap.find(vertexIndex2) != vertexMap.end()){
        mesh::Edge& e = filteredMesh.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
        edgeMap[edge.getID()] = &e;
      }
    }

    //add all triangles formed by contributing edges
    if(_dimensions==3){
      foreach (mesh::Triangle& triangle, seed.triangles() ){
        int edgeIndex1 = triangle.edge(0).getID();
        int edgeIndex2 = triangle.edge(1).getID();
        int edgeIndex3 = triangle.edge(2).getID();
        if(edgeMap.find(edgeIndex1) != edgeMap.end() &&
           edgeMap.find(edgeIndex2) != edgeMap.end() &&
           edgeMap.find(edgeIndex3) != edgeMap.end() ){
          filteredMesh.createTriangle(*edgeMap[edgeIndex1],*edgeMap[edgeIndex2],*edgeMap[edgeIndex3]);
        }
      }
    }

    seed.clear();
    clearBoundingMappings();
    seed.addMesh(filteredMesh);

    preciceDebug("Filtered mesh. #vertices: " << seed.vertices().size()
         <<", #edges: " << seed.edges().size()
         <<", #triangles: " << seed.triangles().size() << ", rank: " << _rank);

    int numberOfVertices = globalVertexIDs.size();
    _masterSlaveCom->send(numberOfVertices,0);
    if(numberOfVertices!=0){
      _masterSlaveCom->send(raw(globalVertexIDs),numberOfVertices,0);
    }

  }
  else{ //master
    assertion(_rank==0);
    assertion(_size>1);
    for(int rankSlave = 1; rankSlave < _size; rankSlave++){
      //slaves have ranks from 0 to size-2
      //TODO better rewrite accept/request connection
      com::CommunicateMesh(_masterSlaveCom).sendMesh ( seed, rankSlave-1 );
      int numberOfVertices = -1;
      _masterSlaveCom->receive(numberOfVertices,rankSlave-1);
      std::vector<int> globalVertexIDs(numberOfVertices,-1);
      if(numberOfVertices!=0){
        _masterSlaveCom->receive(raw(globalVertexIDs),numberOfVertices,rankSlave-1);
      }
      seed.getVertexDistribution()[rankSlave] = globalVertexIDs;
    }
  }
}

void CommunicatedGeometry:: setBoundingFromMapping(
    mapping::PtrMapping mapping)
{
  _boundingFromMapping = mapping;
}

void CommunicatedGeometry:: setBoundingToMapping(
    mapping::PtrMapping mapping)
{
  _boundingToMapping = mapping;
}

void CommunicatedGeometry:: computeBoundingMappings()
{
  if(_boundingFromMapping.use_count()>0){
    _boundingFromMapping->computeMapping();
  }
  if(_boundingToMapping.use_count()>0){
    _boundingToMapping->computeMapping();
  }
}

void CommunicatedGeometry:: clearBoundingMappings()
{
  if(_boundingFromMapping.use_count()>0){
    _boundingFromMapping->clear();
  }
  if(_boundingToMapping.use_count()>0){
    _boundingToMapping->clear();
  }
}

bool CommunicatedGeometry:: doesVertexContribute(
    int vertexID)
{
  //works as easy as this since only read-consistent and write-conservative are allowed
  assertion(_boundingFromMapping.use_count()>0||_boundingToMapping.use_count()>0);
  bool exit = false;
  if(_boundingFromMapping.use_count()>0){
    exit = exit || _boundingFromMapping->doesVertexContribute(vertexID);
  }
  if(_boundingToMapping.use_count()>0){
    exit = exit || _boundingToMapping->doesVertexContribute(vertexID);
  }
  return exit;
}

}} // namespace precice, geometry
