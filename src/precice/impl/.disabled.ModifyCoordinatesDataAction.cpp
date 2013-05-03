// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ModifyCoordinatesDataAction.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "geometry/Geometry.hpp"
#include "spacetree/Spacetree.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/WrappedVector.h"

namespace precice {
namespace impl {

tarch::logging::Log ModifyCoordinatesDataAction::
  _log ( "precice::impl::ModifyCoordinatesDataAction" );


ModifyCoordinatesDataAction:: ModifyCoordinatesDataAction
(
  TimingConstants timing,
  int             dataID,
  MeshContext &   meshContext,
  ModeConstants   mode )
:
  AbstractDataAction ( timing ),
  _data ( meshContext.mesh->data(dataID) ),
  _meshContext ( meshContext ),
  _mode ( mode )
{
  preciceCheck ( _data->getType() == mesh::Data::TYPE_VECTOR,
                 "ModifyCoordinatesMeshAction()",
                 "ModifyCoordinatesMeshAction can only be used for "
                 << "vector type data!"  );
}

//bool ModifyCoordinatesDataAction:: doesUseMeshContext
//(
//   MeshContext * meshContext )
//{
//  foreach ( const mesh::PtrData & data, meshContext->mesh->data() ) {
//    if ( data->getID() == _data->getID() ) {
//      return true;
//    }
//  }
//  return false;
//}

void ModifyCoordinatesDataAction:: performAction
(
  double dt,
  double computedPartFullDt,
  double fullDt )
{
  preciceTrace ( "performAction()" );
  using tarch::la::slice;
  mesh::Data::Values & values = _data->values();
  mesh::PtrMesh mesh = _meshContext.mesh;
  Vector data;
  const int dim = utils::Def::DIM;
  if ( _mode == ADD_TO_COORDINATES_MODE ) {
    preciceDebug ( "Adding data to coordinates" );
    foreach ( mesh::Vertex & vertex, mesh->vertices() ) {
      assign(data) = slice<dim>(values, vertex.getID()*dim);
//      data = tarch::la::dview ( values, vertex.getID() );
      vertex.setCoords ( vertex.getCoords() + data );
    }
  }
  else if ( _mode == SUBTRACT_FROM_COORDINATES_MODE ) {
    preciceDebug ( "Subtracting data from coordinates" );
    foreach ( mesh::Vertex & vertex, mesh->vertices() ) {
      assign(data) = slice<dim>(values, vertex.getID()*dim);
//      data = tarch::la::dview ( values, vertex.getID() );
      vertex.setCoords ( vertex.getCoords() - data );
    }
  }
  else {
    preciceError ( "performAction()", "Unknown mode type!" );
  }
  mesh->computeState ();

  // Reset spacetree, if one is used
  spacetree::PtrSpacetree tree = _meshContext.spacetree;
  preciceCheck ( tree.use_count() == 0, "performAction()",
    "Spacetrees are not yet supported in data actions!");
  if ( tree.use_count() > 0 ) {
    tree->clear ();
    tree->initialize ( *mesh );
  }
}


}} // namespace precice, impl
