#include "cplscheme/CouplingData.hpp"

#include <utility>

#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::cplscheme {

CouplingData::CouplingData(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          requiresInitialization,
    int           extrapolationOrder)
    : requiresInitialization(requiresInitialization),
      _data(std::move(data)),
      _mesh(std::move(mesh)),
      _extrapolation(extrapolationOrder)
{
  PRECICE_ASSERT(_data != nullptr);
  _previousIteration = Eigen::VectorXd::Zero(getSize());
  PRECICE_ASSERT(_mesh != nullptr);
  PRECICE_ASSERT(_mesh.use_count() > 0);
}

int CouplingData::getDimensions() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->getDimensions();
}

int CouplingData::getSize() const
{
  PRECICE_ASSERT(_data != nullptr);
  return values().size();
}

Eigen::VectorXd &CouplingData::values()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->values();
}

const Eigen::VectorXd &CouplingData::values() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->values();
}

Eigen::MatrixXd &CouplingData::gradients()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->gradients();
}

const Eigen::MatrixXd &CouplingData::gradients() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->gradients();
}

bool CouplingData::hasGradient() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->hasGradient();
}

int CouplingData::meshDimensions() const
{
  return _mesh->getDimensions();
}

void CouplingData::storeIteration()
{
  _previousIteration = this->values();
  if (this->hasGradient()) {
    PRECICE_ASSERT(this->gradients().size() > 0);
    _previousIterationGradients = this->gradients();
  }
}

const Eigen::VectorXd CouplingData::previousIteration() const
{
  return _previousIteration;
}

const Eigen::MatrixXd &CouplingData::previousIterationGradients() const
{
  return _previousIterationGradients;
}

int CouplingData::getPreviousIterationSize() const
{
  return previousIteration().size();
}

int CouplingData::getMeshID()
{
  return _mesh->getID();
}

int CouplingData::getDataID()
{
  return _data->getID();
}

std::string CouplingData::getDataName()
{
  return _data->getName();
}

std::vector<int> CouplingData::getVertexOffsets()
{
  return _mesh->getVertexOffsets();
}

void CouplingData::initializeExtrapolation()
{
  _extrapolation.initialize(getSize());
  storeIteration();
}

void CouplingData::moveToNextWindow()
{
  _extrapolation.moveToNextWindow();
  values() = _extrapolation.getInitialGuess();
}

void CouplingData::storeExtrapolationData()
{
  _extrapolation.store(values());
}

} // namespace precice::cplscheme
