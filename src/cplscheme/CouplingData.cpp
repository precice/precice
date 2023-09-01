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
      _previousIteration(_data->getDimensions(), Eigen::VectorXd::Zero(getSize())),
      _exchangeSubsteps(exchangeSubsteps)
{
  PRECICE_ASSERT(_data != nullptr);
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
  this->sample()     = stamples.back().sample;
  _previousIteration = this->sample();
}

const Eigen::VectorXd CouplingData::previousIteration() const
{
  return _previousIteration.values;
}

const Eigen::MatrixXd &CouplingData::previousIterationGradients() const
{
  return _previousIteration.gradients;
}

int CouplingData::getPreviousIterationSize() const
{
  return _previousIteration.values.size();
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

std::string CouplingData::getMeshName()
{
  return _mesh->getName();
}

std::vector<int> CouplingData::getVertexOffsets()
{
  return _mesh->getVertexOffsets();
}

void CouplingData::moveToNextWindow()
{
  _data->moveToNextWindow();
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
