#include "cplscheme/GlobalCouplingData.hpp"

#include <utility>

// #include "mesh/GlobalData.hpp"
#include "mesh/Data.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::cplscheme {

GlobalCouplingData::GlobalCouplingData(
    mesh::PtrData data,
    bool          requiresInitialization,
    int           extrapolationOrder)
    : requiresInitialization(requiresInitialization),
      _data(std::move(data)),
      _extrapolation(extrapolationOrder)
{
  PRECICE_ASSERT(_data != nullptr);
  _previousIteration = Eigen::VectorXd::Zero(getSize());
}

int GlobalCouplingData::getDimensions() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->getDimensions();
}

int GlobalCouplingData::getSize() const
{
  PRECICE_ASSERT(_data != nullptr);
  return values().size();
}

Eigen::VectorXd &GlobalCouplingData::values()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->values();
}

const Eigen::VectorXd &GlobalCouplingData::values() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->values();
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
  return _data->getID();
}

std::string GlobalCouplingData::getDataName()
{
  return _data->getName();
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
