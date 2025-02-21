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
    bool          exchangeSubsteps,
    Direction     direction)
    : requiresInitialization(requiresInitialization),
      _mesh(std::move(mesh)),
      _data(std::move(data)),
      _previousTimeStepsStorage(),
      _exchangeSubsteps(exchangeSubsteps),
      _direction(direction)
{
  PRECICE_ASSERT(_data != nullptr);
  _previousTimeStepsStorage = _data->timeStepsStorage();

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
  // @ŧodo this correct implementation breaks a ton of tests that don't define vertices of a test mesh
  //return _mesh->nVertices() * getDimensions();
  return sample().values.size();
}

int CouplingData::nVertices() const
{
  // @ŧodo this correct implementation breaks a ton of tests that don't define vertices of a test mesh
  //return _mesh->nVertices();
  return sample().values.size() / getDimensions();
}

const Eigen::VectorXd &CouplingData::values() const
{
  return sample().values;
}

const Eigen::MatrixXd &CouplingData::gradients() const
{
  return sample().gradients;
}

int CouplingData::gradientsRows() const
{
  return sample().gradients.rows();
}

int CouplingData::gradientsCols() const
{
  return sample().gradients.cols();
}

const time::Sample &CouplingData::sample() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->timeStepsStorage().last().sample;
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

time::SampleResult CouplingData::getPreviousValuesAtTime(double relativeDt)
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
  this->_sample() = sample; // @todo at some point we should not need this anymore, when mapping, acceleration ... directly work on _timeStepsStorage
  _data->setSampleAtTime(time, sample);
}

void CouplingData::initializeWithZeroAtTime(double time)
{
  this->_sample().setZero(); // @todo at some point we should not need this anymore, when mapping, acceleration ... directly work on _timeStepsStorage
  _data->setSampleAtTime(time, this->sample());
}

void CouplingData::emplaceSampleAtTime(double time)
{
  _data->emplaceSampleAtTime(time);
}

void CouplingData::emplaceSampleAtTime(double time, std::initializer_list<double> values)
{
  this->_sample() = time::Sample(_data->getDimensions(), Eigen::Map<const Eigen::VectorXd>(values.begin(), values.size())); // @todo at some point we should not need this anymore, when mapping, acceleration ... directly work on _timeStepsStorage
  _data->emplaceSampleAtTime(time, values);
}

void CouplingData::emplaceSampleAtTime(double time, std::initializer_list<double> values, std::initializer_list<double> gradients)
{
  this->_sample() = time::Sample(_data->getDimensions(), Eigen::Map<const Eigen::VectorXd>(values.begin(), values.size()), Eigen::Map<const Eigen::MatrixXd>(gradients.begin(), nVertices(), _data->getDimensions() * _data->getSpatialDimensions())); // @todo at some point we should not need this anymore, when mapping, acceleration ... directly work on _timeStepsStorage
  _data->emplaceSampleAtTime(time, values, gradients);
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

void CouplingData::reinitialize()
{
  // TODO port this to subcyling

  // The mesh was reinitialized and new written data will be added later in advance().
  // Meaning all samples are based on a different mesh.
  // Without remapping, the best we can do is setting them to zero samples.
  // We keep the timestamps not to break convergence measures, accelerations, and actions
  auto zero = time::Sample(_data->getDimensions(), _mesh->nVertices());
  zero.setZero();

  _data->timeStepsStorage().setAllSamples(zero);
  _previousTimeStepsStorage.setAllSamples(zero);
}

void CouplingData::storeIteration()
{
  const auto &stamples = this->stamples();
  PRECICE_ASSERT(stamples.size() > 0);
  this->_sample()           = stamples.back().sample;
  _previousTimeStepsStorage = _data->timeStepsStorage();
}

const Eigen::VectorXd &CouplingData::previousIteration() const
{
  PRECICE_ASSERT(!_previousTimeStepsStorage.stamples().empty());
  return _previousTimeStepsStorage.stamples().back().sample.values;
}

const Eigen::MatrixXd &CouplingData::previousIterationGradients() const
{
  PRECICE_ASSERT(!_previousTimeStepsStorage.stamples().empty());
  return _previousTimeStepsStorage.stamples().back().sample.gradients;
}

int CouplingData::getPreviousIterationSize() const
{
  PRECICE_ASSERT(!_previousTimeStepsStorage.stamples().empty());
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

std::string CouplingData::getDataName() const
{
  return _data->getName();
}

std::string CouplingData::getMeshName() const
{
  return _mesh->getName();
}

std::vector<int> CouplingData::getVertexOffsets()
{
  return _mesh->getVertexOffsets();
}

CouplingData::Direction CouplingData::getDirection() const
{
  return _direction;
}

void CouplingData::moveToNextWindow()
{
  if (_direction == Direction::Receive) {
    //_data->moveToNextWindow();
    // _previousTimeStepsStorage = _data->timeStepsStorage();
  }
  _data->moveToNextWindow();
  _previousTimeStepsStorage = _data->timeStepsStorage();
}

time::Sample &CouplingData::_sample()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->sample();
}

bool CouplingData::exchangeSubsteps() const
{
  return _exchangeSubsteps;
}
} // namespace precice::cplscheme
