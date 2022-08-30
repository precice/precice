#include "cplscheme/CouplingData.hpp"

#include <utility>

#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace cplscheme {

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
  _previousIteration = Eigen::VectorXd::Zero(_data->values().size());
  PRECICE_ASSERT(_mesh != nullptr);
  PRECICE_ASSERT(_mesh.use_count() > 0);
}

int CouplingData::getDimensions() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->getDimensions();
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
}

const Eigen::VectorXd CouplingData::previousIteration() const
{
  return _previousIteration;
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
  _extrapolation.initialize(values().size());
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

int CouplingData::getNumberOfStoredTimeSteps()
{
  return _timeStepsStorage.size();
}

Eigen::VectorXd CouplingData::getStoredTimesAscending()
{
  // create std::vector with all keys
  std::vector<double> keys;
  for (auto timeStep : _timeStepsStorage) {
    keys.push_back(timeStep.first);
  }

  // sort vector
  std::sort(keys.begin(), keys.end());

  // copy data into Eigen::VectorXd to return
  auto times = Eigen::VectorXd(keys.size());
  for (int i = 0; i < keys.size(); i++) {
    times[i] = keys[i];
  }
  return times;
}

void CouplingData::clearTimeStepsStorage()
{
  _timeStepsStorage.clear();
}

double CouplingData::maxStoredDt()
{
  double maxDt = -1;
  for (auto timeStep : _timeStepsStorage) {
    if (timeStep.first > maxDt) {
      maxDt = timeStep.first;
    }
  }
  return maxDt;
}

void CouplingData::storeDataAtTime(double relativeDt)
{
  PRECICE_ASSERT(relativeDt > 0.0);
  //PRECICE_ASSERT(relativeDt > maxStoredDt());  // generally a nice security check, but currently we have to override some data after acceleration was performed.
  PRECICE_ASSERT(relativeDt <= 1.0);
  _timeStepsStorage[relativeDt] = _data->values();
}

Eigen::VectorXd CouplingData::getDataAtTime(double relativeDt)
{
  PRECICE_ASSERT(relativeDt > 0.0);
  PRECICE_ASSERT(relativeDt <= 1.0);
  PRECICE_ASSERT(_timeStepsStorage.count(relativeDt) > 0);
  return _timeStepsStorage[relativeDt];
}

} // namespace cplscheme
} // namespace precice
