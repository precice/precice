// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_SPACETREE_IMPL_ENVIRONMENT_HPP_
#define PRECICE_SPACETREE_IMPL_ENVIRONMENT_HPP_

#include "spacetree/Spacetree.hpp"
#include "tarch/Assertions.h"
#include "tarch/la/DynamicVector.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace spacetree {
namespace impl {

class Environment
{
public:

  Environment ( int cellSize, int neighborCellSize )
  :
    _neighborCellIndices(cellSize),
    _neighborSideIndices(cellSize),
    _position(Spacetree::positionUndefined()),
    _neighborCellPositions(neighborCellSize)
  {}

  Environment ( const Environment& toCopy );

  const tarch::la::DynamicVector<int>& getNeighborCellIndices ( int cellIndex ) const
  {
    assertion2((cellIndex >= 0) && (cellIndex < _neighborCellIndices.size()),
               cellIndex, _neighborCellIndices.size());
    return _neighborCellIndices[cellIndex];
  }

  const tarch::la::DynamicVector<int>& getNeighborSideIndices ( int cellIndex ) const
  {
    assertion2((cellIndex >= 0) && (cellIndex < _neighborSideIndices.size()),
               cellIndex, _neighborSideIndices.size());
    return _neighborSideIndices[cellIndex];
  }

  const tarch::la::DynamicVector<int>& getNeighborCellPositions() const
  {
    return _neighborCellPositions;
  }

  int getPosition() const
  {
    return _position;
  }

  void setNeighborCellIndices (
    int                                  cellIndex,
    const tarch::la::DynamicVector<int>& indices )
  {
    assertion2((cellIndex >= 0) && (cellIndex < _neighborCellIndices.size()),
               cellIndex, _neighborCellIndices.size());
    if (_neighborCellIndices[cellIndex].size() == 0){
      _neighborCellIndices[cellIndex].append(indices);
    }
    else {
      _neighborCellIndices[cellIndex] = indices;
    }
  }

  void setNeighborSideIndices (
    int                                  cellIndex,
    const tarch::la::DynamicVector<int>& indices )
  {
    assertion2((cellIndex >= 0) && (cellIndex < _neighborSideIndices.size()),
               cellIndex, _neighborSideIndices.size());
    if (_neighborSideIndices[cellIndex].size() == 0){
      _neighborSideIndices[cellIndex].append(indices);
    }
    else {
      _neighborSideIndices[cellIndex] = indices;
    }
  }

  void setNeighborCellPosition ( int sideIndex, int position )
  {
    assertion2((sideIndex >= 0) && (sideIndex < _neighborCellPositions.size()),
               sideIndex, _neighborCellPositions.size());
    _neighborCellPositions[sideIndex] = position;
  }

  void setAllNeighborCellPositions ( int position )
  {
    assign(_neighborCellPositions) = position;
  }

  void computePosition();



private:

  static tarch::logging::Log _log;

  tarch::la::DynamicVector<tarch::la::DynamicVector<int> > _neighborCellIndices;

  tarch::la::DynamicVector<tarch::la::DynamicVector<int> > _neighborSideIndices;

  int _position;

  tarch::la::DynamicVector<int> _neighborCellPositions;
};

}}} // namespace precice, spacetree, impl

#endif /* PRECICE_SPACETREE_IMPL_ENVIRONMENT_HPP_ */
