#include <boost/range.hpp>

#include "cplscheme/CouplingScheme.hpp"
#include "math/Bspline.hpp"
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

void Storage::setSampleAtTime(double time, const Sample &sample)
{
  // The spline has to be recomputed, since the underlying data has changed
  _bspline.reset();

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
  } else {
    // Overriding sample
    existingSample->sample = sample;
  }
}

void Storage::setAllSamples(const Sample &sample)
{
  for (auto &stample : _stampleStorage) {
    stample.sample = sample;
  }
}

void Storage::setInterpolationDegree(int interpolationDegree)
{
  PRECICE_ASSERT(interpolationDegree >= Time::MIN_WAVEFORM_DEGREE);
  _degree = interpolationDegree;

  // The spline has to be recomputed, since the underlying data has changed
  _bspline.reset();
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
  _stampleStorage.erase(_stampleStorage.begin(), --_stampleStorage.end());
  PRECICE_ASSERT(nextWindowStart == _stampleStorage.front().timestamp);

  // The spline has to be recomputed, since the underlying data has changed
  _bspline.reset();
}

void Storage::trim()
{
  PRECICE_ASSERT(!_stampleStorage.empty(), "Storage does not contain any data!");
  const double thisWindowStart = _stampleStorage.front().timestamp;
  _stampleStorage.erase(++_stampleStorage.begin(), _stampleStorage.end());
  PRECICE_ASSERT(_stampleStorage.size() == 1);
  PRECICE_ASSERT(thisWindowStart == _stampleStorage.front().timestamp);

  // The spline has to be recomputed, since the underlying data has changed
  _bspline.reset();
}

void Storage::clear()
{
  _stampleStorage.clear();
  PRECICE_ASSERT(_stampleStorage.size() == 0);

  // The spline has to be recomputed, since the underlying data has changed
  _bspline.reset();
}

void Storage::clearExceptLast()
{
  if (_stampleStorage.empty()) {
    return;
  }
  _stampleStorage.erase(_stampleStorage.begin(), --_stampleStorage.end());

  // The spline has to be recomputed, since the underlying data has changed
  _bspline.reset();
}

void Storage::trimBefore(double time)
{
  auto beforeTime = [time](const auto &s) { return math::smaller(s.timestamp, time); };
  _stampleStorage.erase(std::remove_if(_stampleStorage.begin(), _stampleStorage.end(), beforeTime), _stampleStorage.end());

  // The spline has to be recomputed, since the underlying data has changed
  _bspline.reset();
}

void Storage::trimAfter(double time)
{
  auto afterTime = [time](const auto &s) { return math::greater(s.timestamp, time); };
  _stampleStorage.erase(std::remove_if(_stampleStorage.begin(), _stampleStorage.end(), afterTime), _stampleStorage.end());

  // The spline has to be recomputed, since the underlying data has changed
  _bspline.reset();
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

bool Storage::empty() const
{
  return _stampleStorage.empty();
}

const time::Stample &Storage::last() const
{
  PRECICE_ASSERT(!_stampleStorage.empty());
  return _stampleStorage[_stampleStorage.size() - 1];
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
  PRECICE_ASSERT(this->nTimes() != 0, "There are no samples available");
  const int usedDegree = computeUsedDegree(_degree, nTimes());

  if (usedDegree == 0) {
    return this->getSampleAtOrAfter(time).values;
  }

  PRECICE_ASSERT(usedDegree >= 1);

  // Find existing samples
  for (const auto &stample : _stampleStorage) {
    if (math::equals(stample.timestamp, time)) {
      return stample.sample.values;
    }
    if (math::greater(stample.timestamp, time)) {
      break;
    }
  }

  //Create a new bspline if _bspline does not already contain a spline
  if (!_bspline.has_value()) {
    auto [times, values] = getTimesAndValues();
    _bspline.emplace(times, values, usedDegree);
  }
  return _bspline.value().interpolateAt(time);
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
  return std::min(requestedDegree, numberOfAvailableSamples - 1);
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
