#include "com/SerializedStamples.hpp"
#include "cplscheme/CouplingData.hpp"
#include "math/differences.hpp"

namespace precice::com::serialize {

SerializedStamples SerializedStamples::serialize(const cplscheme::PtrCouplingData data)
{
  SerializedStamples result;

  result._timeSteps = 2; // @todo use all available time steps for subcycling

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

void SerializedStamples::deserializeInto(Eigen::VectorXd timeStamps, const cplscheme::PtrCouplingData data)
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
  const auto &atBeginn = data->timeStepsStorage().stamples().front();
  int         timeId   = 0;
  PRECICE_ASSERT(atBeginn.timestamp == time::Storage::WINDOW_START, atBeginn.timestamp);
  auto sliceBeginn = atBeginn.sample.values;
  for (int valueId = 0; valueId < data->getSize(); valueId++) {
    _values(valueId * _timeSteps) = sliceBeginn(valueId);
  }

  const auto &atEnd = data->timeStepsStorage().stamples().back();
  timeId            = 1;
  if (data->timeStepsStorage().nTimes() == 1) { // during exchangeInitialData where only data at WINDOW_START is provided. Use identical data atBeginn and atEnd.
    PRECICE_ASSERT(math::equals(atEnd.timestamp, time::Storage::WINDOW_START), atEnd.timestamp);
  } else {
    PRECICE_ASSERT(math::equals(atEnd.timestamp, time::Storage::WINDOW_END), atEnd.timestamp);
  }
  auto sliceEnd = atEnd.sample.values;
  for (int valueId = 0; valueId < data->getSize(); valueId++) {
    _values(valueId * _timeSteps + timeId) = sliceEnd(valueId);
  }
}

void SerializedStamples::serializeGradients(const cplscheme::PtrCouplingData data)
{
  const auto &atBeginn = data->timeStepsStorage().stamples().front();
  int         timeId   = 0;
  PRECICE_ASSERT(atBeginn.timestamp == time::Storage::WINDOW_START, atBeginn.timestamp);
  const auto sliceBeginn = Eigen::VectorXd::Map(atBeginn.sample.gradients.data(), atBeginn.sample.gradients.rows() * atBeginn.sample.gradients.cols());
  PRECICE_ASSERT(data->getSize() * data->meshDimensions() * data->getDimensions() == sliceBeginn.size(), data->getSize(), data->meshDimensions(), data->getDimensions(), sliceBeginn.size());
  for (int valueId = 0; valueId < sliceBeginn.size(); valueId++) {
    _gradients(valueId * _timeSteps) = sliceBeginn(valueId);
  }

  const auto &atEnd = data->timeStepsStorage().stamples().back();
  timeId            = 1;
  if (data->timeStepsStorage().nTimes() == 1) { // during exchangeInitialData where only data at WINDOW_START is provided. Use identical data atBeginn and atEnd.
    PRECICE_ASSERT(math::equals(atEnd.timestamp, time::Storage::WINDOW_START), atEnd.timestamp);
  } else {
    PRECICE_ASSERT(math::equals(atEnd.timestamp, time::Storage::WINDOW_END), atEnd.timestamp);
  }
  const auto sliceEnd = Eigen::VectorXd::Map(atEnd.sample.gradients.data(), atEnd.sample.gradients.rows() * atEnd.sample.gradients.cols());
  PRECICE_ASSERT(data->getSize() * data->meshDimensions() * data->getDimensions() == sliceEnd.size(), data->getSize(), data->meshDimensions(), data->getDimensions(), sliceEnd.size());
  for (int valueId = 0; valueId < data->getSize() * data->meshDimensions(); valueId++) {
    _gradients(valueId * _timeSteps + timeId) = sliceEnd(valueId);
  }
}

void SerializedStamples::deserialize(const Eigen::VectorXd timeStamps, cplscheme::PtrCouplingData data) const
{
  PRECICE_ASSERT(timeStamps.size() * data->getSize() == _values.size(), timeStamps.size() * data->getSize(), _values.size());

  data->timeStepsStorage().trim();

  for (int timeId = 0; timeId < timeStamps.size(); timeId++) {
    const double time = timeStamps(timeId);
    PRECICE_ASSERT(math::greaterEquals(time, time::Storage::WINDOW_START) && math::greaterEquals(time::Storage::WINDOW_END, time)); // time < 0 or time > 1 is not allowed.

    Eigen::VectorXd slice(data->getSize());
    for (int valueId = 0; valueId < slice.size(); valueId++) {
      slice(valueId) = _values(valueId * timeStamps.size() + timeId);
    }

    if (not data->hasGradient()) {
      data->setSampleAtTime(time, time::Sample{slice});
    } else {
      PRECICE_ASSERT(data->hasGradient());

      Eigen::MatrixXd gradientSlice(data->sample().gradients.rows(), data->sample().gradients.cols());
      auto            gradientView = Eigen::VectorXd::Map(gradientSlice.data(), gradientSlice.rows() * gradientSlice.cols());
      for (int gradientId = 0; gradientId < gradientView.size(); gradientId++) {
        gradientView(gradientId) = _gradients(gradientId * timeStamps.size() + timeId);
      }
      data->setSampleAtTime(time, time::Sample{slice, gradientSlice});
    }
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
