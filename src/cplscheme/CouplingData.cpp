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
      _mesh(std::move(mesh)), _extrapolationOrder(extrapolationOrder),
      _timeStepsStorageCurrent(extrapolationOrder), _timeStepsStoragePrevious(extrapolationOrder)
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

  _timeStepsStoragePrevious.clear(false);

  // Need a better way of copying the contents of waveform iterations
  int  nValues        = getSize();
  int  nTimeSteps     = _timeStepsStorageCurrent.nTimes();
  auto timesAndValues = _timeStepsStorageCurrent.getTimesAndValues();
  auto serializedData = timesAndValues.second;
  auto timesAscending = timesAndValues.first;

  for (int timeId = 0; timeId < timesAscending.size(); timeId++) {
    auto slice = Eigen::VectorXd(getSize());
    for (int valueId = 0; valueId < slice.size(); valueId++) {
      slice(valueId) = serializedData(valueId * timesAscending.size() + timeId);
    }
    auto time = timesAscending(timeId);
    PRECICE_ASSERT(math::greaterEquals(time, time::Storage::WINDOW_START) && math::greaterEquals(time::Storage::WINDOW_END, time)); // time < 0 or time > 1 is not allowed.
    _timeStepsStoragePrevious.setValuesAtTime(time, slice);
  }

  this->_previousIteration = this->getValuesAtTime(_timeStepsStorageCurrent.maxStoredNormalizedDt()); // store data in values to make this non-breaking.

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

void CouplingData::initializeStorage(Eigen::VectorXd data)
{
  _timeStepsStorageCurrent.clear(false); // only required for MultiCouplingScheme, if participant that is not the controller has send data (with data) which is also non-initialized receive data. See DataBC in Integration/Serial/MultiCoupling/MultiCouplingThreeSolvers3.
  storeValuesAtTime(time::Storage::WINDOW_START, data);
  storeValuesAtTime(time::Storage::WINDOW_END, data);
}

Eigen::VectorXd CouplingData::getStoredTimesAscending()
{
  return _timeStepsStorageCurrent.getTimes();
}

void CouplingData::clearTimeStepsStorage()
{
  _timeStepsStorageCurrent.clear();
}

void CouplingData::moveTimeStepsStorage()
{
  _timeStepsStorageCurrent.move();
  values() = _timeStepsStorageCurrent.getValuesAtTime(time::Storage::WINDOW_END); // @todo Better do this just before returning to SolverInterfaceImpl. Compare to BaseCouplingScheme::receiveData
}

void CouplingData::storeValuesAtTime(double relativeDt, Eigen::VectorXd data, bool mustOverwriteExisting)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  _timeStepsStorageCurrent.setValuesAtTime(relativeDt, data, mustOverwriteExisting);
}

Eigen::VectorXd CouplingData::getValuesAtTime(double relativeDt)
{
  return _timeStepsStorageCurrent.getValuesAtTime(relativeDt);
}

Eigen::VectorXd CouplingData::getPreviousValuesAtTime(double relativeDt)
{

  return _timeStepsStoragePrevious.getValuesAtOrAfter(relativeDt);
}

Eigen::VectorXd CouplingData::getSerialized()
{
  int  nValues        = getSize();
  int  nTimeSteps     = _timeStepsStorageCurrent.nTimes();
  auto serializedData = Eigen::VectorXd(nTimeSteps * nValues);
  auto timesAndValues = _timeStepsStorageCurrent.getTimesAndValues();
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

  _timeStepsStorageCurrent.clear(false);

  for (int timeId = 0; timeId < timesAscending.size(); timeId++) {
    auto slice = Eigen::VectorXd(getSize());
    for (int valueId = 0; valueId < slice.size(); valueId++) {
      slice(valueId) = serializedData(valueId * timesAscending.size() + timeId);
    }
    auto time = timesAscending(timeId);
    PRECICE_ASSERT(math::greaterEquals(time, time::Storage::WINDOW_START) && math::greaterEquals(time::Storage::WINDOW_END, time)); // time < 0 or time > 1 is not allowed.
    this->storeValuesAtTime(time, slice);
  }
  this->values() = this->getValuesAtTime(_timeStepsStorageCurrent.maxStoredNormalizedDt()); // store data in values to make this non-breaking.
}

} // namespace precice::cplscheme
