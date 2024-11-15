#include "com/SerializedStamples.hpp"
#include "cplscheme/CouplingData.hpp"
#include "math/differences.hpp"

namespace precice::com::serialize {

SerializedStamples SerializedStamples::serialize(const cplscheme::PtrCouplingData data)
{
  SerializedStamples result;

  result._timeSteps = data->timeStepsStorage().nTimes(); // @todo use all available time steps for subcycling
  result.allocate(data);
  result.serializeValues(data);
  if (data->hasGradient()) {
    result.serializeGradients(data);
  }
  return result;
}

SerializedStamples SerializedStamples::empty(Eigen::VectorXd timeStamps, const cplscheme::PtrCouplingData data)
{
  SerializedStamples result;

  result._timeSteps = timeStamps.size();

  result.allocate(data);

  return result;
}

void SerializedStamples::deserializeInto(const Eigen::VectorXd &timeStamps, const cplscheme::PtrCouplingData data)
{
  PRECICE_ASSERT(_timeSteps == timeStamps.size());

  deserialize(timeStamps, data);
}

void SerializedStamples::allocate(const cplscheme::PtrCouplingData data)
{
  _values = Eigen::VectorXd(_timeSteps * data->getSize());

  if (data->hasGradient()) {
    _gradients = Eigen::VectorXd(_timeSteps * data->getSize() * data->meshDimensions() * data->getDimensions());
  }
}

void SerializedStamples::serializeValues(const cplscheme::PtrCouplingData data)
{
  const int nValues = data->sample().values.size();
  int       timeId  = 0;
  for (const auto &stample : data->stamples()) {
    const Eigen::VectorXd &slice = stample.sample.values;
    for (int valueId = 0; valueId < nValues; valueId++) {
      _values(valueId * _timeSteps + timeId) = slice(valueId);
    }
    timeId++;
  }
}

void SerializedStamples::serializeGradients(const cplscheme::PtrCouplingData data)
{
  const int nValues = data->sample().gradients.size();
  int       timeId  = 0;
  for (const auto &stample : data->stamples()) {
    const Eigen::VectorXd &slice = Eigen::VectorXd::Map(stample.sample.gradients.data(), stample.sample.gradients.rows() * stample.sample.gradients.cols());
    PRECICE_ASSERT(nValues == slice.size());
    for (int valueId = 0; valueId < slice.size(); valueId++) {
      _gradients(valueId * _timeSteps + timeId) = slice(valueId);
    }
    timeId++;
  }
}

void SerializedStamples::deserialize(const Eigen::VectorXd timeStamps, cplscheme::PtrCouplingData data) const
{
  PRECICE_ASSERT(timeStamps.size() * data->getSize() == _values.size(), timeStamps.size() * data->getSize(), _values.size());

  data->timeStepsStorage().clear(); // @todo needs optimization. Don't need to communicate and serialize / deserialize data at beginning of window, because it is already there.

  const auto dataDims = data->getDimensions();

  for (int timeId = 0; timeId < timeStamps.size(); timeId++) {
    const double time = timeStamps(timeId);

    Eigen::VectorXd slice(data->getSize());
    for (int valueId = 0; valueId < slice.size(); valueId++) {
      slice(valueId) = _values(valueId * timeStamps.size() + timeId);
    }

    if (!data->hasGradient()) {
      data->setSampleAtTime(time, time::Sample{dataDims, std::move(slice)});
      continue;
    }

    Eigen::MatrixXd gradientSlice(data->sample().gradients.rows(), data->sample().gradients.cols());
    auto            gradientView = Eigen::VectorXd::Map(gradientSlice.data(), gradientSlice.rows() * gradientSlice.cols());
    for (int gradientId = 0; gradientId < gradientView.size(); gradientId++) {
      gradientView(gradientId) = _gradients(gradientId * timeStamps.size() + timeId);
    }
    data->setSampleAtTime(time, time::Sample{dataDims, std::move(slice), std::move(gradientSlice)});
  }
}

const Eigen::VectorXd &SerializedStamples::values() const
{
  return _values;
}

Eigen::VectorXd &SerializedStamples::values()
{
  return _values;
}

const Eigen::VectorXd &SerializedStamples::gradients() const
{
  return _gradients;
}

Eigen::VectorXd &SerializedStamples::gradients()
{
  return _gradients;
}

int SerializedStamples::nTimeSteps() const
{
  return _timeSteps;
}

} // namespace precice::com::serialize
