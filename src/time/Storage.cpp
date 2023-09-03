#include <boost/range.hpp>

#include "cplscheme/CouplingScheme.hpp"
#include "math/differences.hpp"
#include "time/Storage.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

Storage::Storage()
    : _stampleStorage{}
{
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
  // @todo reactivate this assertion!
  // PRECICE_ASSERT(nTimes() >= 2, "Calling Storage::move() is only allowed, if there is a sample at the beginning and at the end. This ensures that this function is only called at the end of the window.", getTimes());

  PRECICE_ASSERT(!_stampleStorage.empty(), "Storage does not contain any data!");
  const double nextWindowStart = _stampleStorage.back().timestamp;
  _stampleStorage.erase(_stampleStorage.begin(), --_stampleStorage.end());
  PRECICE_ASSERT(nextWindowStart == _stampleStorage.front().timestamp);
  PRECICE_DEBUG("times after: {}", getTimes());
}

void Storage::trim()
{
  PRECICE_DEBUG("Storage::trim()");
  PRECICE_DEBUG("times before: {}", getTimes());
  PRECICE_ASSERT(!_stampleStorage.empty(), "Storage does not contain any data!");
  const double thisWindowStart = _stampleStorage.front().timestamp;
  _stampleStorage.erase(++_stampleStorage.begin(), _stampleStorage.end());
  PRECICE_ASSERT(_stampleStorage.size() == 1);
  PRECICE_ASSERT(thisWindowStart == _stampleStorage.front().timestamp);
  PRECICE_DEBUG("times after: {}", getTimes());
}

void Storage::clear()
{
  _stampleStorage.clear();
  PRECICE_ASSERT(_stampleStorage.size() == 0);
}

Eigen::VectorXd Storage::getValuesAtOrAfter(double before) const
{
  PRECICE_DEBUG("getValuesAtOrAfter({})", before);
  PRECICE_DEBUG("available times: {}", getTimes());
  if (nTimes() == 1) {
    return _stampleStorage.front().sample.values; // @todo in this case the name getValuesAtOrAfter does not fit, because _stampleStorage.front().sample.values is returned for any time before.
  } else {
    auto stample = std::find_if(_stampleStorage.begin(), _stampleStorage.end(), [&before](const auto &s) { return math::greaterEquals(s.timestamp, before); });
    PRECICE_ASSERT(stample != _stampleStorage.end(), "no values found!");
    return stample->sample.values;
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

time::Sample Storage::getSampleAtBeginning()
{
  return _stampleStorage.front().sample;
}

time::Sample Storage::getSampleAtEnd()
{
  return _stampleStorage.back().sample;
}

} // namespace precice::time
