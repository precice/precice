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

void Storage::initialize(time::Sample sample)
{
  _currentWindowStart = 0;
  _stampleStorage.emplace_back(Stample{_currentWindowStart, sample});
}

void Storage::setSampleAtTime(double time, Sample sample)
{
  PRECICE_ASSERT(not sample.values.hasNaN());
  PRECICE_ASSERT(math::smallerEquals(_currentWindowStart, time), "Setting sample outside of valid range!", _currentWindowStart, time);
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
  PRECICE_ASSERT(!_stampleStorage.empty(), "Storage does not contain any data!");
  PRECICE_ASSERT(math::equals(_stampleStorage.front().timestamp, _currentWindowStart), _stampleStorage.front().timestamp, _currentWindowStart);
  _currentWindowStart = _stampleStorage.back().timestamp;
  _stampleStorage.erase(_stampleStorage.begin(), --_stampleStorage.end());
  PRECICE_ASSERT(_currentWindowStart == _stampleStorage.front().timestamp);
}

void Storage::trim()
{
  PRECICE_ASSERT(!_stampleStorage.empty(), "Storage does not contain any data!");
  PRECICE_ASSERT(math::equals(_stampleStorage.front().timestamp, _currentWindowStart), _stampleStorage.front().timestamp, _currentWindowStart);
  _stampleStorage.erase(++_stampleStorage.begin(), _stampleStorage.end());
}

Eigen::VectorXd Storage::getValuesAtOrAfter(double before) const
{
  if (nTimes() == 1) {
    return _stampleStorage.front().sample.values;
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
