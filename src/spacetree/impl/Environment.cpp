#include "spacetree/impl/Environment.hpp"

namespace precice {
namespace spacetree {
namespace impl {

logging::Logger Environment:: _log("precice::spacetree::impl::Environment");

Environment::Environment ( int cellSize, int neighborCellSize )
  :
    _neighborCellIndices(cellSize),
    _neighborSideIndices(cellSize),
    _position(Spacetree::positionUndefined()),
    _neighborCellPositions(neighborCellSize)
  {}


Environment:: Environment
(
  const Environment& toCopy )
:
  _neighborCellIndices(toCopy._neighborCellIndices),
  _neighborSideIndices(toCopy._neighborSideIndices),
  _position(toCopy._position),
  _neighborCellPositions(toCopy._neighborCellPositions)
{}

const Eigen::VectorXi& Environment::getNeighborCellIndices ( int cellIndex ) const
{
  assertion((cellIndex >= 0) && (cellIndex < _neighborCellIndices.size()),
            cellIndex, _neighborCellIndices.size());
  return _neighborCellIndices[cellIndex];
}

const Eigen::VectorXi& Environment::getNeighborSideIndices ( int cellIndex ) const
{
  assertion((cellIndex >= 0) && (cellIndex < _neighborSideIndices.size()),
            cellIndex, _neighborSideIndices.size());
  return _neighborSideIndices[cellIndex];
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
