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

time::Storage &CouplingData::timeStepsStorage()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->timeStepsStorage();
}

const Eigen::VectorXd &CouplingData::values() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->values();
}

Eigen::MatrixXd &CouplingData::gradientValues()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->gradientValues();
}

const Eigen::MatrixXd &CouplingData::gradientValues() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->gradientValues();
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
  const auto stamples = this->timeStepsStorage().getStamples();
  // PRECICE_ASSERT(stamples.size() > 0);  //@todo preferable, but cannot use this, because of some invalid configs in tests (e.g. tests/serial/AitkenAcceleration.xml)
  if (stamples.size() > 0) {
    this->values() = stamples.back().sample.values;
    if (this->hasGradient()) {
      this->gradientValues() = stamples.back().sample.gradient;
    }
  }

  _previousIteration = this->values();
  if (this->hasGradient()) {
    PRECICE_ASSERT(this->gradientValues().size() > 0);
    _previousIterationGradients = this->gradientValues();
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

  if (this->hasGradient()) {
    this->timeStepsStorage().setSampleAtTime(time::Storage::WINDOW_END, time::Sample{this->values(), this->gradientValues()});
  }
  this->timeStepsStorage().setSampleAtTime(time::Storage::WINDOW_END, time::Sample{this->values()});
}

void CouplingData::storeExtrapolationData()
{
  const auto stamples = this->timeStepsStorage().getStamples();
  // PRECICE_ASSERT(stamples.size() > 0);  //@todo preferable, but cannot use this, because of some invalid configs in tests (e.g. tests/serial/AitkenAcceleration.xml)
  if (stamples.size() > 0) {
    this->values() = stamples.back().sample.values;
  }

  _extrapolation.store(values());
}

} // namespace precice::cplscheme
