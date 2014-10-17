// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CommunicatedGeometry.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
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
      //TODO dann testen 1 nastin, 3 solidz
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

    _boundingFromMapping->computeMapping();
    _boundingToMapping->computeMapping();

    mesh::Mesh boundingMesh("BoundingMesh", _dimensions, seed.isFlipNormals());

    foreach ( const mesh::Vertex& vertex, seed.vertices() ){
      boundingMesh.createVertex(vertex.getCoords());
    }
    //TODO auch edges erzeugen etc.
    seed.clear();

    preciceDebug("Before: Bounding Mesh vertices: " << boundingMesh.vertices().size()
        << ", seed vertices: " << seed.vertices().size());

    std::vector<int> globalVertexIDs;
    int verticesSize = boundingMesh.vertices().size();
    for(int i=0; i < verticesSize; i++){
      if(_boundingFromMapping->doesVertexContribute(i)||_boundingToMapping->doesVertexContribute(i)){
        seed.createVertex(boundingMesh.vertices()[i].getCoords());
        globalVertexIDs.push_back(i);
        //TODO hier wirds komplizierter, nur edges zwischen 2 contributing hinzufügen
        //TODO auch doesVertexContr für projection mapping implementieren
      }
    }
    _boundingFromMapping->clear();
    _boundingToMapping->clear();

    preciceDebug("After: Bounding Mesh vertices: " << boundingMesh.vertices().size()
        << ", seed vertices: " << seed.vertices().size());

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

}} // namespace precice, geometry
