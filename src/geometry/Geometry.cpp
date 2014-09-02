// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Geometry.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/Globals.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
//#include "utils/NumericalCompare.hpp"
#include "boost/foreach.hpp"

namespace precice {
namespace geometry {

tarch::logging::Log Geometry:: _log ( "precice::geometry::Geometry" );

Geometry:: Geometry
(
  const utils::DynVector& offset )
:
  _vertexDistribution(),
  _isDistributed(false),
  _isProvided(false),
  _offset (offset)
{}

void Geometry:: create
(
  mesh::Mesh& seed )
{
  preciceTrace1 ( "create()", seed.getName() );
  assertion2 ( seed.getDimensions() == _offset.size(), seed.getDimensions(),
               _offset.size() );
  specializedCreate ( seed );
  utils::DynVector zero ( seed.getDimensions(), 0.0 );
  if ( not tarch::la::equals(getOffset(), zero) ) {
    utils::DynVector temp ( seed.getDimensions() );
    foreach ( mesh::Vertex& vertex, seed.vertices() ) {
      temp = _offset;
      temp += vertex.getCoords();
      vertex.setCoords ( temp );
    }
  }
  seed.computeState();
  allocateDataValues(seed);
}

void Geometry:: allocateDataValues ( mesh::Mesh & mesh )
{
  mesh.allocateDataValues ();
}

void Geometry:: collectDistribution(
    mesh::Mesh& seed,
    com::PtrCommunication   masterSlaveCom,
    const int               rank,
    const int               size)
{
  preciceTrace1 ( "collectDistribution()", rank );
  _isDistributed = true;

  if(rank>0){ //slave
    com::CommunicateMesh(masterSlaveCom).sendMesh ( seed, 0 );
  }
  else{ //master
    assertion(rank==0);
    assertion(size>1);
    _vertexDistribution[0]=seed.vertices().size();
    mesh::Mesh slaveMesh("SlaveMesh", seed.getDimensions(), seed.isFlipNormals());
    mesh::Mesh& rSlaveMesh = slaveMesh;

    for(int rankSlave = 1; rankSlave < size; rankSlave++){
      //slaves have ranks from 0 to size-2
      //TODO better rewrite accept/request connection
      rSlaveMesh.clear();
      com::CommunicateMesh(masterSlaveCom).receiveMesh ( rSlaveMesh, rankSlave-1);
      _vertexDistribution[rankSlave]=rSlaveMesh.vertices().size();
      int dim = rSlaveMesh.getDimensions();
      utils::DynVector coord(dim);
      foreach ( const mesh::Vertex& vertex, rSlaveMesh.vertices() ){
          coord = vertex.getCoords();
          seed.createVertex(coord);
      }
    }
  }
}




}} // namespace precice, geometry
