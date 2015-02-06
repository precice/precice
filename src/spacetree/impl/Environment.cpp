// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "spacetree/impl/Environment.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace spacetree {
namespace impl {

tarch::logging::Log Environment:: _log("precice::spacetree::impl::Environment");

Environment:: Environment
(
  const Environment& toCopy )
:
  _neighborCellIndices(toCopy._neighborCellIndices.size()),
  _neighborSideIndices(toCopy._neighborSideIndices.size()),
  _position(toCopy._position),
  _neighborCellPositions(toCopy._neighborCellPositions)
{
  //preciceTrace ( "Environment(Environment)" );
  for ( int i=0; i < _neighborCellIndices.size(); i++ ){
    assertion1 ( toCopy._neighborCellIndices[i].size() > 0, i );
    _neighborCellIndices[i].append(toCopy._neighborCellIndices[i]);
  }
  for ( int i=0; i < _neighborSideIndices.size(); i++ ){
    assertion1 ( toCopy._neighborSideIndices[i].size() > 0, i );
    _neighborSideIndices[i].append(toCopy._neighborSideIndices[i]);
  }
}

void Environment:: computePosition()
{
  for ( int i=0; i < (int)_neighborCellPositions.size(); i++ ){
    int pos = _neighborCellPositions[i];
    if (pos == Spacetree::positionUndefined()){
      continue;
    }
    else if (pos == Spacetree::positionOnGeometry()){
      _position = pos;
    }
    else {
      _position = pos;
      return;
    }
  }
  preciceCheck( _position == Spacetree::positionOnGeometry(),
                "computePosition()", "Could compute position!" );
}

}}} // namespace precice, spacetree, impl
