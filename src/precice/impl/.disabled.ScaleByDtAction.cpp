// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ScaleByDtAction.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/PropertyContainer.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace impl {

tarch::logging::Log ScaleByDtAction:: _log ( "precice::impl::ScaleByDtAction" );

ScaleByDtAction:: ScaleByDtAction
(
  TimingConstants       timing,
  int                   sourceDataID,
  int                   targetDataID,
  const mesh::PtrMesh & mesh,
  Scaling               scaling )
:
  AbstractDataAction ( timing ),
  _sourceData ( mesh->data(sourceDataID) ),
  _targetData ( mesh->data (targetDataID) ),
  _mesh ( mesh ),
  _scaling ( scaling )
{
  assertion2 ( _sourceData->getType() == _targetData->getType(),
               _sourceData->getTypeName(), _targetData->getTypeName() );
}

void ScaleByDtAction:: performAction
(
  double dt,
  double computedPartFullDt,
  double fullDt )
{
  preciceTrace ( "performAction()" );
  mesh::Data::Values & sourceValues = _sourceData->values();
  mesh::Data::Values & targetValues = _targetData->values();
  assertion2 ( sourceValues.size() == targetValues.size(),
               sourceValues.size(), targetValues.size() );
  if ( _scaling == SCALING_BY_COMPUTED_TIMESTEP ){
    double scaling = dt / fullDt;
    for ( int i=0; i < targetValues.size(); i++ ){
      targetValues[i] = sourceValues[i] * scaling;
    }
  }
  else {
    assertion1 ( _scaling == SCALING_BY_COMPUTED_TIMESTEP_PART, _scaling );
    double scaling = computedPartFullDt / fullDt;
    for ( int i=0; i < targetValues.size(); i++ ){
      targetValues[i] = sourceValues[i] * scaling;
    }
  }
}

}} // namspace precice, impl
