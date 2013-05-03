// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ScaleByAreaAction.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace impl {

tarch::logging::Log ScaleByAreaAction:: _log ( "precice::impl::ScaleByAreaAction" );

ScaleByAreaAction:: ScaleByAreaAction
(
  TimingConstants       timing,
  int                   targetDataID,
  const mesh::PtrMesh & mesh,
  Scaling               scaling )
:
  AbstractDataAction ( timing ),
  _targetData ( mesh->data (targetDataID) ),
  _mesh ( mesh ),
  _scaling ( scaling )
{}

void ScaleByAreaAction:: performAction
(
  double dt,
  double computedPartFullDt,
  double fullDt )
{
  preciceTrace ( "performAction()" );
  preciceCheck ( utils::Def::DIM == 2, "performAction()",
                 "Not implemented for dim != 2!" );
  mesh::Data::Values & targetValues = _targetData->values();
  mesh::Data::Values areas ( _mesh->vertices().size(), 0.0 );
  foreach ( mesh::Edge & edge, _mesh->edges() ) {
    areas[edge.vertex(0).getID()] += edge.getEnclosingRadius();
    areas[edge.vertex(1).getID()] += edge.getEnclosingRadius();
  }
  int dimensions = _targetData->getValueDimension ();
  assertion ( targetValues.size() / dimensions == areas.size() );
  if ( _scaling == SCALING_DIVIDE_BY_AREA ){
    for ( int i=0; i < areas.size(); i++ ) {
      for ( int dim=0; dim < dimensions; dim++ ) {
        int valueIndex = i*dimensions + dim;
        targetValues[valueIndex] /= areas[i];
      }
    }
  }
  else if ( _scaling == SCALING_MULTIPLY_BY_AREA ){
    for ( int i=0; i < areas.size(); i++ ) {
      for ( int dim=0; dim < dimensions; dim++ ) {
        int valueIndex = i*dimensions + dim;
        targetValues[valueIndex] *= areas[i];
      }
    }
  }
}

}} // namspace precice, impl
