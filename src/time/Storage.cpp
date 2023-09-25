#include <boost/range.hpp>

#include "cplscheme/CouplingScheme.hpp"
#include "math/bspline.hpp"
#include "math/differences.hpp"
#include "time/Storage.hpp"
#include "time/Time.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

Storage::Storage()
    : _stampleStorage{}, _degree(0)
{
}

Storage &Storage::operator=(const Storage &other)
{
  this->clear();
  this->_degree = other.getInterpolationDegree();
  for (const auto &stample : other.stamples()) {
    this->setSampleAtTime(stample.timestamp, stample.sample);
  }
  return *this;
}

void Storage::setSampleAtTime(double time, Sample sample)
{
  if (_stampleStorage.empty()) {
    _stampleStorage.emplace_back(Stample{time, sample});
    return;
  }

  const double currentWindowStart = _stampleStorage.front().timestamp;

  PRECICE_ASSERT(not sample.values.hasNaN());
  PRECICE_ASSERT(math::smallerEquals(currentWindowStart, time), "Setting sample outside of valid range!", currentWindowStart, time);
  // check if key "time" exists.
  auto existingSample = std::find_if(_stampleStorage.begin(), _stampleStorage.end(), [&time](const auto &s) { return math::equals(s.timestamp, time); });
  if (existingSample == _stampleStorage.end()) { // key does not exist yet
    PRECICE_ASSERT(math::smaller(maxStoredTime(), time), maxStoredTime(), time, "Trying to write sample with a time that is too small. Please use clear(), if you want to write new samples to the storage.");
    _stampleStorage.emplace_back(Stample{time, sample});
  } else { // overwrite sample at "time"
    for (auto &stample : _stampleStorage) {
      if (math::equals(stample.timestamp, time)) {
        stample.sample = sample;
        return;
      }
    }
    PRECICE_ASSERT(false, "unreachable!");
  }
}

void Storage::setInterpolationDegree(int interpolationDegree)
{
  PRECICE_ASSERT(Time::MIN_WAVEFORM_DEGREE <= _degree && _degree <= Time::MAX_WAVEFORM_DEGREE);
  _degree = interpolationDegree;
}

int Storage::getInterpolationDegree() const
{
  return _degree;
}

double Storage::maxStoredTime() const
{
  if (_stampleStorage.size() == 0) {
    return -1; // invalid return
  } else {
    return _stampleStorage.back().timestamp;
  }
}

int Storage::nTimes() const
{
  return _stampleStorage.size();
}

int Storage::nDofs() const
{
  PRECICE_ASSERT(_stampleStorage.size() > 0);
  return _stampleStorage[0].sample.values.size();
}

void Storage::move()
{
  PRECICE_ASSERT(nTimes() >= 2, "Calling Storage::move() is only allowed, if there is a sample at the beginning and at the end. This ensures that this function is only called at the end of the window.", getTimes());
  PRECICE_ASSERT(!_stampleStorage.empty(), "Storage does not contain any data!");
  const double nextWindowStart = _stampleStorage.back().timestamp;
  clearBefore(nextWindowStart);
  PRECICE_ASSERT(nextWindowStart == _stampleStorage.front().timestamp);
}

void Storage::trim()
{
  PRECICE_ASSERT(!_stampleStorage.empty(), "Storage does not contain any data!");
  const double thisWindowStart = _stampleStorage.front().timestamp;
  clearAfter(thisWindowStart);
  PRECICE_ASSERT(_stampleStorage.size() == 1);
  PRECICE_ASSERT(thisWindowStart == _stampleStorage.front().timestamp);
}

void Storage::clearBefore(double time)
{
  PRECICE_DEBUG("Storage::clearBefore({})", time);
  PRECICE_DEBUG("_stampleStorage.size() before is {}", _stampleStorage.size());
  auto stample = std::find_if(_stampleStorage.begin(), _stampleStorage.end(), [&time](const auto &s) { return math::greaterEquals(s.timestamp, time); });
  PRECICE_ASSERT(stample != _stampleStorage.end(), "no values found!", getTimes(), time);
  _stampleStorage.erase(_stampleStorage.begin(), stample);
  PRECICE_DEBUG("_stampleStorage.size() after is {}", _stampleStorage.size());
}

void Storage::clearAfter(double time)
{
  auto stample = std::find_if(_stampleStorage.begin(), _stampleStorage.end(), [&time](const auto &s) { return math::greater(s.timestamp, time); });
  if (stample != _stampleStorage.end()) { // @todo remove this safeguard?
    PRECICE_ASSERT(stample != _stampleStorage.end(), "no values found!");
    _stampleStorage.erase(stample, _stampleStorage.end());
  } else {
    // PRECICE_ASSERT(math::equals(_stampleStorage.back().timestamp, time));
  }
}

void Storage::clear()
{
  _stampleStorage.clear();
  PRECICE_ASSERT(_stampleStorage.size() == 0);
}

Sample Storage::getSampleAtOrAfter(double before) const
{
  PRECICE_TRACE(before);
  if (nTimes() == 1) {
    return _stampleStorage.front().sample; // @todo in this case the name getSampleAtOrAfter does not fit, because _stampleStorage.front().sample is returned for any time before.
  } else {
    auto stample = std::find_if(_stampleStorage.begin(), _stampleStorage.end(), [&before](const auto &s) { return math::greaterEquals(s.timestamp, before); });
    PRECICE_ASSERT(stample != _stampleStorage.end(), "no values found!");
    return stample->sample;
  }
}

Eigen::VectorXd Storage::getTimes() const
{
  auto times = Eigen::VectorXd(nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i] = _stampleStorage[i].timestamp;
  }
  return times;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> Storage::getTimesAndValues() const
{
  auto times  = Eigen::VectorXd(nTimes());
  auto values = Eigen::MatrixXd(nDofs(), nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i]      = _stampleStorage[i].timestamp;
    values.col(i) = _stampleStorage[i].sample.values;
  }
  return std::make_pair(times, values);
}

Eigen::VectorXd Storage::sample(double time) const
{
  const int usedDegree = computeUsedDegree(_degree, nTimes());

  if (usedDegree == 0) {
    return this->getSampleAtOrAfter(time).values;
  }

  PRECICE_ASSERT(usedDegree >= 1);

  auto data = getTimesAndValues();

  return math::bspline::interpolateAt(data.first, data.second, usedDegree, time);
}

Eigen::MatrixXd Storage::sampleGradients(double time) const
{
  const int usedDegree = computeUsedDegree(_degree, nTimes());

  if (usedDegree == 0) {
    return this->getSampleAtOrAfter(time).gradients;
  }

  PRECICE_WARN("You specified interpolation degree of {}, but only degree 0 is supported for gradient interpolation", usedDegree); // @todo implement this like for sampleAt
  return this->getSampleAtOrAfter(time).gradients;
}

int Storage::computeUsedDegree(int requestedDegree, int numberOfAvailableSamples) const
{
  int usedDegree = -1;
  PRECICE_ASSERT(requestedDegree <= 3);
  if (requestedDegree == 0 || numberOfAvailableSamples < 2) {
    usedDegree = 0;
  } else if (requestedDegree == 1 || numberOfAvailableSamples < 3) {
    usedDegree = 1;
  } else if (requestedDegree == 2 || numberOfAvailableSamples < 4) {
    usedDegree = 2;
  } else if (requestedDegree == 3 || numberOfAvailableSamples < 5) {
    usedDegree = 3;
  } else {
    PRECICE_ASSERT(false); // not supported
  }
  return usedDegree;
}

time::Sample Storage::getSampleAtBeginning()
{
  return _stampleStorage.front().sample;
}

time::Sample Storage::getSampleAtEnd()
{
  return _stampleStorage.back().sample;
}

} // namespace precice::time
