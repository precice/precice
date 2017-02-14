#pragma once

#include "utils/assertion.hpp"
#include "logging/Logger.hpp"
#include <Eigen/Core>

namespace precice {
namespace spacetree {
namespace impl {

class Environment
{
public:

  Environment ( int cellSize, int neighborCellSize );

  Environment ( const Environment& toCopy );

  const Eigen::VectorXi& getNeighborCellIndices ( size_t cellIndex ) const;

  const Eigen::VectorXi& getNeighborSideIndices ( size_t cellIndex ) const;
  
  const Eigen::VectorXi& getNeighborCellPositions() const
  {
    return _neighborCellPositions;
  }

  int getPosition() const
  {
    return _position;
  }

  void setNeighborCellIndices (
    size_t cellIndex,
    const Eigen::VectorXi& indices )
  {
    assertion(cellIndex < _neighborCellIndices.size(), cellIndex, _neighborCellIndices.size());
    _neighborCellIndices[cellIndex] = indices;
  }

  void setNeighborSideIndices (
    size_t cellIndex,
    const Eigen::VectorXi& indices )
  {
    assertion(cellIndex < _neighborSideIndices.size(),
              cellIndex, _neighborSideIndices.size());
    _neighborSideIndices[cellIndex] = indices;
  }

  void setNeighborCellPosition ( size_t sideIndex, int position )
  {
    assertion(sideIndex < static_cast<size_t>(_neighborCellPositions.size()),
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

