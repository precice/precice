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
    int           extrapolationOrder)
    : requiresInitialization(requiresInitialization),
      _exchangeSubsteps(exchangeSubsteps),
      _data(std::move(data)),
      _mesh(std::move(mesh))
{
  PRECICE_ASSERT(_data != nullptr);
  /// Lazy allocation of _previousIteration.gradient: only used in case the corresponding data has gradients
  _previousIteration = time::Sample{Eigen::VectorXd::Zero(getSize())};
  timeStepsStorage().setExtrapolationOrder(extrapolationOrder);

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
  this->sample() = sample; // @todo at some point we should not need this anymore, when mapping, acceleration ... directly work on _timeStepsStorage
  timeStepsStorage().setSampleAtTime(time, sample);
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

std::vector<int> CouplingData::getVertexOffsets()
{
  return _mesh->getVertexOffsets();
}

void CouplingData::moveToNextWindow()
{
  if (this->timeStepsStorage().stamples().size() > 0) {
    this->timeStepsStorage().move();
    const auto &atEnd = this->timeStepsStorage().stamples().back();
    PRECICE_ASSERT(math::equals(atEnd.timestamp, time::Storage::WINDOW_END));
    _data->sample() = atEnd.sample;
  }
}

Eigen::VectorXd CouplingData::getSerializedValues()
{
  const int nValues    = getSize();
  int       nTimeSteps = timeStepsStorage().nTimes();

  if (nTimeSteps == 1) { // special treatment during initialization
    PRECICE_ASSERT(timeStepsStorage().stamples().front().timestamp == time::Storage::WINDOW_START);
    nTimeSteps = 2;
    Eigen::VectorXd serializedData(nTimeSteps * nValues);
    int             timeId = 0;

    for (int i = 0; i < nTimeSteps; i++) { // put sample at WINDOW_START twice into serialized data
      const auto &           stample = timeStepsStorage().stamples().front();
      const Eigen::VectorXd &slice   = stample.sample.values;
      for (int valueId = 0; valueId < nValues; valueId++) {
        serializedData(valueId * nTimeSteps + timeId) = slice(valueId);
      }
      timeId++;
    }

    return serializedData;
  } else { // normal treatment where every sample is serialized
    Eigen::VectorXd serializedData(nTimeSteps * nValues);
    int             timeId = 0;

    for (const auto &stample : timeStepsStorage().stamples()) {
      const Eigen::VectorXd &slice = stample.sample.values;
      for (int valueId = 0; valueId < nValues; valueId++) {
        serializedData(valueId * nTimeSteps + timeId) = slice(valueId);
      }
      timeId++;
    }
    return serializedData;
  }
}

Eigen::VectorXd CouplingData::getSerializedGradients()
{
  const int nValues    = sample().gradients.size();
  int       nTimeSteps = timeStepsStorage().nTimes();

  if (nTimeSteps == 1) { // special treatment during initialization
    PRECICE_ASSERT(timeStepsStorage().stamples().front().timestamp == time::Storage::WINDOW_START);
    nTimeSteps = 2;
    Eigen::VectorXd serializedGradientsData(nTimeSteps * nValues);
    int             timeId = 0;

    for (int i = 0; i < nTimeSteps; i++) { // put sample at WINDOW_START twice into serialized data
      const auto &           stample = timeStepsStorage().stamples().front();
      const Eigen::VectorXd &slice   = Eigen::VectorXd::Map(stample.sample.gradients.data(), stample.sample.gradients.rows() * stample.sample.gradients.cols());
      PRECICE_ASSERT(nValues == slice.size());
      for (int valueId = 0; valueId < slice.size(); valueId++) {
        serializedGradientsData(valueId * nTimeSteps + timeId) = slice(valueId);
      }
      timeId++;
    }

    return serializedGradientsData;

  } else { // normal treatment where every sample is serialized
    Eigen::VectorXd serializedGradientsData(nTimeSteps * nValues);
    int             timeId = 0;

    for (const auto &stample : timeStepsStorage().stamples()) {
      const Eigen::VectorXd &slice = Eigen::VectorXd::Map(stample.sample.gradients.data(), stample.sample.gradients.rows() * stample.sample.gradients.cols());
      PRECICE_ASSERT(nValues == slice.size());
      for (int valueId = 0; valueId < slice.size(); valueId++) {
        serializedGradientsData(valueId * nTimeSteps + timeId) = slice(valueId);
      }
      timeId++;
    }

    return serializedGradientsData;
  }
}

void CouplingData::storeFromSerialized(Eigen::VectorXd timesAscending, Eigen::VectorXd serializedValues)
{
  PRECICE_ASSERT(timesAscending.size() * getSize() == serializedValues.size());

  timeStepsStorage().trim();

  for (int timeId = 0; timeId < timesAscending.size(); timeId++) {
    Eigen::VectorXd slice(getSize());
    for (int valueId = 0; valueId < slice.size(); valueId++) {
      slice(valueId) = serializedValues(valueId * timesAscending.size() + timeId);
    }
    const double time = timesAscending(timeId);
    PRECICE_ASSERT(math::greaterEquals(time, time::Storage::WINDOW_START) && math::greaterEquals(time::Storage::WINDOW_END, time)); // time < 0 or time > 1 is not allowed.
    this->setSampleAtTime(time, time::Sample{slice});
  }
}

void CouplingData::storeFromSerialized(Eigen::VectorXd timesAscending, Eigen::VectorXd serializedValues, Eigen::MatrixXd serializedGradients)
{
  PRECICE_ASSERT(timesAscending.size() * getSize() == serializedValues.size());

  timeStepsStorage().trim();

  for (int timeId = 0; timeId < timesAscending.size(); timeId++) {
    Eigen::VectorXd slice(getSize());
    for (int valueId = 0; valueId < slice.size(); valueId++) {
      slice(valueId) = serializedValues(valueId * timesAscending.size() + timeId);
    }

    Eigen::MatrixXd gradientSlice(sample().gradients.rows(), sample().gradients.cols());
    auto            gradientView = Eigen::VectorXd::Map(gradientSlice.data(), gradientSlice.rows() * gradientSlice.cols());
    for (int gradientId = 0; gradientId < gradientView.size(); gradientId++) {
      gradientView(gradientId) = serializedGradients(gradientId * timesAscending.size() + timeId);
    }
    const double time = timesAscending(timeId);
    PRECICE_ASSERT(math::greaterEquals(time, time::Storage::WINDOW_START) && math::greaterEquals(time::Storage::WINDOW_END, time)); // time < 0 or time > 1 is not allowed.
    this->setSampleAtTime(time, time::Sample{slice, gradientSlice});
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
