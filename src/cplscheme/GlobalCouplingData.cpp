#include "cplscheme/GlobalCouplingData.hpp"

#include <utility>

#include "mesh/GlobalData.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::cplscheme {

GlobalCouplingData::GlobalCouplingData(
    mesh::PtrGlobalData globalData,
    bool                requiresInitialization,
    int                 extrapolationOrder)
    : requiresInitialization(requiresInitialization),
      _globalData(std::move(globalData)),
      _extrapolation(extrapolationOrder)
{
  PRECICE_ASSERT(_globalData != nullptr);
  _previousIteration = Eigen::VectorXd::Zero(getSize());
}

int GlobalCouplingData::getDimensions() const
{
  PRECICE_ASSERT(_globalData != nullptr);
  return _globalData->getDimensions();
}

int GlobalCouplingData::getSize() const
{
  PRECICE_ASSERT(_globalData != nullptr);
  return values().size();
}

Eigen::VectorXd &GlobalCouplingData::values()
{
  PRECICE_ASSERT(_globalData != nullptr);
  return _globalData->values();
}

const Eigen::VectorXd &GlobalCouplingData::values() const
{
  PRECICE_ASSERT(_globalData != nullptr);
  return _globalData->values();
}

void GlobalCouplingData::storeIteration()
{
  _previousIteration = this->values();
}

const Eigen::VectorXd GlobalCouplingData::previousIteration() const
{
  return _previousIteration;
}

int GlobalCouplingData::getPreviousIterationSize() const
{
  return previousIteration().size();
}

int GlobalCouplingData::getDataID()
{
  return _globalData->getID();
}

std::string GlobalCouplingData::getDataName()
{
  return _globalData->getName();
}

void GlobalCouplingData::initializeExtrapolation()
{
  _extrapolation.initialize(getSize());
  storeIteration();
}

void GlobalCouplingData::moveToNextWindow()
{
  _extrapolation.moveToNextWindow();
  values() = _extrapolation.getInitialGuess();
}

void GlobalCouplingData::storeExtrapolationData()
{
  _extrapolation.store(values());
}

} // namespace precice::cplscheme
