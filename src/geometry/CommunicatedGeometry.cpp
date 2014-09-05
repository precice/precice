// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CommunicatedGeometry.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "mesh/Mesh.hpp"
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
  _dimensions(dimensions)
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
      seed.setDistributionType(mesh::Mesh::EXACT);
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
      seed.setDistributionType(mesh::Mesh::ALL);
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
      foreach ( const mesh::Vertex& vertex, rSlaveMesh.vertices() ){
          coord = vertex.getCoords();
          seed.createVertex(coord);
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

  if(_rank>0){ //slave
    com::CommunicateMesh(_masterSlaveCom).receiveMesh ( seed, 0);
  }
  else{ //master
    assertion(_rank==0);
    assertion(_size>1);
    for(int rankSlave = 1; rankSlave < _size; rankSlave++){
      //slaves have ranks from 0 to size-2
      //TODO better rewrite accept/request connection
      //at the moment complete mesh is sent to every slave
      com::CommunicateMesh(_masterSlaveCom).sendMesh ( seed, rankSlave-1 );
      int numberOfVertices = 0;
      foreach ( const mesh::Vertex& vertex, seed.vertices() ){
        seed.getVertexDistribution()[rankSlave].push_back(numberOfVertices);
        numberOfVertices++;
      }
    }
  }
}

}} // namespace precice, geometry
