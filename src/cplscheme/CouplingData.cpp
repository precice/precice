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
      _mesh(std::move(mesh))
{
  // @todo store extrapolation order in Waveform?
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

void CouplingData::moveToNextWindow()
{
  _timeStepsStorage.move();
}

Eigen::VectorXd CouplingData::getStoredTimesAscending()
{
  return _timeStepsStorage.getTimes();
}

void CouplingData::clearTimeStepsStorage(bool keepZero)
{
  _timeStepsStorage.clear(keepZero);
}

void CouplingData::moveTimeStepsStorage()
{
  _timeStepsStorage.move();
}

void CouplingData::storeDataAtTime(Eigen::VectorXd data, double relativeDt)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(relativeDt, _timeStepsStorage.maxStoredNormalizedDt()), relativeDt, _timeStepsStorage.maxStoredNormalizedDt());
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  _timeStepsStorage.setValueAtTime(relativeDt, data);
}

void CouplingData::overrideDataAtEndWindowTime(Eigen::VectorXd data)
{
  _timeStepsStorage.overrideDataAtEndWindowTime(data);
}

Eigen::VectorXd CouplingData::getDataAtTime(double relativeDt)
{
  return _timeStepsStorage.getValueAtTime(relativeDt);
}

Eigen::VectorXd CouplingData::getSerialized()
{
  int  nValues        = getSize();
  int  nTimeSteps     = _timeStepsStorage.nTimes();
  auto serializedData = Eigen::VectorXd(nTimeSteps * nValues);
  auto timesAndValues = _timeStepsStorage.getTimesAndValues();
  auto values         = timesAndValues.second;

  for (int timeId = 0; timeId < nTimeSteps; timeId++) {
    auto slice = values.col(timeId);
    for (int valueId = 0; valueId < nValues; valueId++) {
      serializedData(valueId * nTimeSteps + timeId) = slice(valueId);
    }
  }
  return serializedData;
}

void CouplingData::storeFromSerialized(Eigen::VectorXd timesAscending, Eigen::VectorXd serializedData)
{
  PRECICE_ASSERT(timesAscending.size() * getSize() == serializedData.size());
  for (int timeId = 0; timeId < timesAscending.size(); timeId++) {
    auto slice = Eigen::VectorXd(getSize());
    for (int valueId = 0; valueId < slice.size(); valueId++) {
      slice(valueId) = serializedData(valueId * timesAscending.size() + timeId);
    }
    auto time = timesAscending(timeId);
    PRECICE_ASSERT(math::greaterEquals(time, time::Storage::WINDOW_START) && math::greaterEquals(time::Storage::WINDOW_END, time)); // time < 0 or time > 1 is not allowed.
    this->storeDataAtTime(slice, time);
  }
  this->values() = this->getDataAtTime(_timeStepsStorage.maxStoredNormalizedDt()); // store data in values to make this non-breaking.
}

} // namespace precice::cplscheme
