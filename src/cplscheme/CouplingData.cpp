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
    bool          exchangeSubsteps)
    : requiresInitialization(requiresInitialization),
      _mesh(std::move(mesh)),
      _data(std::move(data)),
      _exchangeSubsteps(exchangeSubsteps),
      _previousTimeStepsStorage()
{
  PRECICE_ASSERT(_data != nullptr);
  _previousTimeStepsStorage.setInterpolationDegree(3); // @todo hard-coded for now, but we need to somehow link this to <read-data waveform-order="ORDER" />
  _previousTimeStepsStorage.setSampleAtTime(time::Storage::WINDOW_START, time::Sample{getDimensions(), Eigen::VectorXd::Zero(getSize())});
  _previousTimeStepsStorage.setSampleAtTime(time::Storage::WINDOW_END, time::Sample{getDimensions(), Eigen::VectorXd::Zero(getSize())});

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
  return sample().values.size();
}

Eigen::VectorXd &CouplingData::values()
{
  return sample().values;
}

const Eigen::VectorXd &CouplingData::values() const
{
  return sample().values;
}

Eigen::MatrixXd &CouplingData::gradients()
{
  return sample().gradients;
}

const Eigen::MatrixXd &CouplingData::gradients() const
{
  return sample().gradients;
}

time::Storage &CouplingData::timeStepsStorage()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->timeStepsStorage();
}

const time::Storage &CouplingData::timeStepsStorage() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->timeStepsStorage();
}

Eigen::VectorXd CouplingData::getPreviousValuesAtTime(double relativeDt)
{
  return _previousTimeStepsStorage.sample(relativeDt);
}

Eigen::MatrixXd CouplingData::getPreviousGradientsAtTime(double relativeDt)
{
  return _previousTimeStepsStorage.sampleGradients(relativeDt);
}

void CouplingData::setSampleAtTime(double time, time::Sample sample)
{
  PRECICE_ASSERT(not sample.values.hasNaN());
  this->sample() = sample; // @todo at some point we should not need this anymore, when mapping, acceleration ... directly work on _timeStepsStorage
  _data->setSampleAtTime(time, sample);
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
  const auto &stamples = this->stamples();
  PRECICE_ASSERT(stamples.size() > 0);
  this->sample() = stamples.back().sample;

  _previousTimeStepsStorage.trim();
  if (stamples.size() == 1) { // special treatment during initialization
    const auto &stample = this->stamples().back();
    PRECICE_ASSERT(math::equals(stample.timestamp, time::Storage::WINDOW_START), "stample.timestamp must be WINDOW_START", stample.timestamp);
    _previousTimeStepsStorage.setSampleAtTime(time::Storage::WINDOW_START, stample.sample);
    _previousTimeStepsStorage.setSampleAtTime(time::Storage::WINDOW_END, stample.sample);
  } else {
    PRECICE_ASSERT(math::equals(this->stamples().back().timestamp, time::Storage::WINDOW_END), "Only allowed to storeIteration, if at window end");
    _previousTimeStepsStorage = _data->timeStepsStorage();
  }
}

const Eigen::VectorXd CouplingData::previousIteration() const
{
  return _previousTimeStepsStorage.stamples().back().sample.values;
}

const Eigen::MatrixXd &CouplingData::previousIterationGradients() const
{
  return _previousTimeStepsStorage.stamples().back().sample.gradients;
}

int CouplingData::getPreviousIterationSize() const
{
  return _previousTimeStepsStorage.stamples().back().sample.values.size();
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

void CouplingData::moveToNextWindow()
{
  _data->moveToNextWindow();
  if (this->stamples().size() > 0) {
    PRECICE_ASSERT(math::equals(this->stamples().back().timestamp, time::Storage::WINDOW_END), "this->stamples() needs to be initialized properly for WINDOW_START and WINDOW_END");
    _previousTimeStepsStorage = _data->timeStepsStorage();
  }
}

time::Sample &CouplingData::sample()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->sample();
}

const time::Sample &CouplingData::sample() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->sample();
}

bool CouplingData::exchangeSubsteps() const
{
  return _exchangeSubsteps;
}

} // namespace precice::cplscheme
