#include <boost/range.hpp>

#include "cplscheme/CouplingScheme.hpp"
#include "math/bspline.hpp"
#include "math/differences.hpp"
#include "time/Storage.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

const double Storage::WINDOW_START = 0.0;

const double Storage::WINDOW_END = 1.0;

Storage::Storage()
    : _stampleStorage{}, _degree(0)
{
}

Storage &Storage::operator=(const Storage &other)
{
  this->_stampleStorage.clear();
  for (const auto &stample : other.stamples()) {
    this->setSampleAtTime(stample.timestamp, stample.sample);
  }
  return *this;
}

void Storage::initialize(time::Sample sample)
{
  _stampleStorage.emplace_back(Stample{WINDOW_START, sample});
  _stampleStorage.emplace_back(Stample{WINDOW_END, sample});
}

void Storage::setSampleAtTime(double time, Sample sample)
{
  PRECICE_ASSERT(not sample.values.hasNaN());
  PRECICE_ASSERT(math::smallerEquals(WINDOW_START, time), "Setting sample outside of valid range!");
  PRECICE_ASSERT(math::smallerEquals(time, WINDOW_END), "Setting sample outside of valid range!");
  // check if key "time" exists.
  auto existingSample = std::find_if(_stampleStorage.begin(), _stampleStorage.end(), [&time](const auto &s) { return math::equals(s.timestamp, time); });
  if (existingSample == _stampleStorage.end()) { // key does not exist yet
    PRECICE_ASSERT(math::smaller(maxStoredNormalizedDt(), time), maxStoredNormalizedDt(), time, "Trying to write sample with a time that is too small. Please use clear(), if you want to write new samples to the storage.");
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
  _degree = interpolationDegree;
}

double Storage::maxStoredNormalizedDt() const
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
  PRECICE_ASSERT(math::equals(_stampleStorage.front().timestamp, time::Storage::WINDOW_START), _stampleStorage.front().timestamp);
  PRECICE_ASSERT(math::equals(_stampleStorage.back().timestamp, time::Storage::WINDOW_END), _stampleStorage.back().timestamp);

  auto sampleAtBeginning = getSampleAtEnd();
  auto sampleAtEnd       = getSampleAtEnd();
  _stampleStorage.clear();
  _stampleStorage.emplace_back(time::Stample{WINDOW_START, std::move(sampleAtBeginning)});
  _stampleStorage.emplace_back(time::Stample{WINDOW_END, std::move(sampleAtEnd)});
}

void Storage::trim()
{
  PRECICE_ASSERT(!_stampleStorage.empty(), "Storage does not contain any data!");
  PRECICE_ASSERT(_stampleStorage.front().timestamp == time::Storage::WINDOW_START);
  _stampleStorage.erase(++_stampleStorage.begin(), _stampleStorage.end());
}

Sample Storage::getSampleAtOrAfter(double before) const
{
  auto stample = std::find_if(_stampleStorage.begin(), _stampleStorage.end(), [&before](const auto &s) { return math::greaterEquals(s.timestamp, before); });
  PRECICE_ASSERT(stample != _stampleStorage.end(), "no values found!");

  return stample->sample;
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

Eigen::VectorXd Storage::sample(double normalizedDt) const
{
  const int usedDegree = computeUsedDegree(_degree, nTimes());

  PRECICE_ASSERT(math::equals(this->maxStoredNormalizedDt(), time::Storage::WINDOW_END), this->maxStoredNormalizedDt()); // sampling is only allowed, if a window is complete.

  if (usedDegree == 0) {
    return this->getSampleAtOrAfter(normalizedDt).values;
  }

  PRECICE_ASSERT(usedDegree >= 1);

  auto data = getTimesAndValues();

  return math::bspline::interpolateAt(data.first, data.second, usedDegree, normalizedDt);
}

Eigen::MatrixXd Storage::sampleGradients(double normalizedDt) const
{
  const int usedDegree = computeUsedDegree(_degree, nTimes());

  PRECICE_ASSERT(math::equals(this->maxStoredNormalizedDt(), time::Storage::WINDOW_END), this->maxStoredNormalizedDt()); // sampling is only allowed, if a window is complete.

  if (usedDegree == 0) {
    return this->getSampleAtOrAfter(normalizedDt).gradients;
  }

  PRECICE_WARN("You specified interpolation degree of {}, but only degree 0 is supported for gradient interpolation", usedDegree); // @todo implement this like for sampleAt
  return this->getSampleAtOrAfter(normalizedDt).gradients;
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
