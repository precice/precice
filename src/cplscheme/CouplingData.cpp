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
  _previousTimeStepsStorage = _data->waveform();

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
  return _mesh->nVertices() * getDimensions();
}

int CouplingData::nVertices() const
{
  return _mesh->nVertices();
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
  const int rows = meshDimensions();
  PRECICE_ASSERT(sample().gradients.rows() == rows, sample().gradients.rows(), rows);
  return rows;
}

int CouplingData::gradientsCols() const
{
  const int cols = getSize();
  PRECICE_ASSERT(sample().gradients.cols() == cols, sample().gradients.cols(), cols);
  return cols;
}

const time::Sample &CouplingData::sample() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->sample();
}

time::Waveform &CouplingData::waveform()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->waveform();
}

const time::Waveform &CouplingData::waveform() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->waveform();
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
  _data->setSampleAtTime(time, sample);
}

void CouplingData::setGlobalSample(const time::Sample &sample)
{
  PRECICE_ASSERT(not sample.values.hasNaN());
  _data->setGlobalSample(sample);
}

void CouplingData::initializeWithZeroAtTime(double time)
{
  if (!hasGradient()) {
    auto zero = time::Sample(getDimensions(), nVertices());
    zero.setZero();
    _data->setSampleAtTime(time, zero);
    return;
  }
  auto zero = time::Sample(getDimensions(), nVertices(), meshDimensions());
  zero.setZero();
  _data->setSampleAtTime(time, zero);
}

void CouplingData::emplaceSampleAtTime(double time)
{
  _data->emplaceSampleAtTime(time);
}

void CouplingData::emplaceSampleAtTime(double time, std::initializer_list<double> values)
{
  _data->emplaceSampleAtTime(time, values);
}

void CouplingData::emplaceSampleAtTime(double time, std::initializer_list<double> values, std::initializer_list<double> gradients)
{
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
  auto zero = time::Sample(getDimensions(), nVertices());
  zero.setZero();

  _data->waveform().setAllSamples(zero);
  _previousTimeStepsStorage.setAllSamples(zero);
}

void CouplingData::storeIteration()
{
  const auto &stamples = this->stamples();
  PRECICE_ASSERT(stamples.size() > 0);
  _previousTimeStepsStorage = _data->waveform();
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
    // _previousTimeStepsStorage = _data->waveform();
  }
  _data->moveToNextWindow();
  _previousTimeStepsStorage = _data->waveform();
}

bool CouplingData::exchangeSubsteps() const
{
  return _exchangeSubsteps;
}
} // namespace precice::cplscheme
