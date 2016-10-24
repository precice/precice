#pragma once

#include "spacetree/Spacetree.hpp"
#include "utils/assertion.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace spacetree {
namespace impl {

class Environment
{
public:

  Environment ( int cellSize, int neighborCellSize );

  Environment ( const Environment& toCopy );

  const Eigen::VectorXi& getNeighborCellIndices ( int cellIndex ) const;

  const Eigen::VectorXi& getNeighborSideIndices ( int cellIndex ) const;
  
  const Eigen::VectorXi& getNeighborCellPositions() const
  {
    return _neighborCellPositions;
  }

  int getPosition() const
  {
    return _position;
  }

  void setNeighborCellIndices (
    int                                  cellIndex,
    const Eigen::VectorXi& indices )
  {
    assertion((cellIndex >= 0) && (cellIndex < _neighborCellIndices.size()),
               cellIndex, _neighborCellIndices.size());
    _neighborCellIndices[cellIndex] = indices;
  }

  void setNeighborSideIndices (
    int                                  cellIndex,
    const Eigen::VectorXi& indices )
  {
    assertion((cellIndex >= 0) && (cellIndex < _neighborSideIndices.size()),
               cellIndex, _neighborSideIndices.size());
    _neighborSideIndices[cellIndex] = indices;
  }

  void setNeighborCellPosition ( int sideIndex, int position )
  {
    assertion((sideIndex >= 0) && (sideIndex < _neighborCellPositions.size()),
               sideIndex, _neighborCellPositions.size());
    _neighborCellPositions[sideIndex] = position;
  }

  void setAllNeighborCellPositions ( int position )
  {
    _neighborCellPositions.setConstant(position);
  }

  void computePosition();



private:

  static logging::Logger _log;

  std::vector<Eigen::VectorXi> _neighborCellIndices;

  std::vector<Eigen::VectorXi> _neighborSideIndices;

  int _position;

  Eigen::VectorXi _neighborCellPositions;
};

}}} // namespace precice, spacetree, impl

