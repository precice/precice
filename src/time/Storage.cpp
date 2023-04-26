#include "time/Storage.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

const double Storage::WINDOW_START = 0.0;

const double Storage::WINDOW_END = 1.0;

Storage::Storage()
    : _sampleStorage{}
{
}

void Storage::initialize(Eigen::VectorXd values)
{
  _sampleStorage.emplace_back(Stample{WINDOW_START, Sample{values}});
  _sampleStorage.emplace_back(Stample{WINDOW_END, Sample{values}});
}

void Storage::setValuesAtTime(double time, Eigen::VectorXd values)
{
  setSampleAtTime(time, Sample{values});
}

void Storage::setSampleAtTime(double time, Sample sample)
{
  PRECICE_ASSERT(math::smallerEquals(WINDOW_START, time), "Setting values outside of valid range!");
  PRECICE_ASSERT(math::smallerEquals(time, WINDOW_END), "Sampling outside of valid range!");
  // check if key "time" exists.
  auto existingSample = std::find_if(_sampleStorage.begin(), _sampleStorage.end(), [&time](const auto &s) { return math::equals(s.timestamp, time); });
  if (existingSample == _sampleStorage.end()) { // key does not exist yet
    PRECICE_ASSERT(math::smaller(maxStoredNormalizedDt(), time), maxStoredNormalizedDt(), time, "Trying to write values with a time that is too small. Please use clear(), if you want to reset the storage.");
    _sampleStorage.emplace_back(Stample{time, sample});
  } else { // overwrite values at "time"
    for (auto &stample : _sampleStorage) {
      if (math::equals(stample.timestamp, time)) {
        stample.sample = sample;
        return;
      }
    }
    PRECICE_ASSERT(false, "unreachable!");
  }
}

double Storage::maxStoredNormalizedDt()
{
  if (_sampleStorage.size() == 0) {
    return -1; // invalid return
  } else {
    return _sampleStorage.back().timestamp;
  }
}

int Storage::nTimes()
{
  return _sampleStorage.size();
}

int Storage::nDofs()
{
  PRECICE_ASSERT(_sampleStorage.size() > 0);
  return _sampleStorage[0].sample.values.size();
}

void Storage::move()
{
  PRECICE_ASSERT(nTimes() > 0);
  auto initialGuess = _sampleStorage.back().sample.values; // use values at end of window as initial guess for next
  _sampleStorage.clear();
  initialize(initialGuess);
}

void Storage::clear()
{
  PRECICE_ASSERT(nTimes() > 0, "Storage does not contain any data!");
  Stample keep = _sampleStorage.front();
  PRECICE_ASSERT(keep.timestamp == time::Storage::WINDOW_START);
  _sampleStorage.clear();
  _sampleStorage.emplace_back(keep);
}

void Storage::clearAll()
{
  _sampleStorage.clear();
}

Eigen::VectorXd Storage::getValuesAtOrAfter(double before)
{
  auto sample = std::find_if(_sampleStorage.begin(), _sampleStorage.end(), [&before](const auto &s) { return math::greaterEquals(s.timestamp, before); });
  PRECICE_ASSERT(sample != _sampleStorage.end(), "no values found!");

  return sample->sample.values;
}

Eigen::VectorXd Storage::getTimes()
{
  auto times = Eigen::VectorXd(nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i] = _sampleStorage[i].timestamp;
  }
  return times;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> Storage::getTimesAndValues()
{
  auto times  = Eigen::VectorXd(nTimes());
  auto values = Eigen::MatrixXd(nDofs(), nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i]      = _sampleStorage[i].timestamp;
    values.col(i) = _sampleStorage[i].sample.values;
  }
  return std::make_pair(times, values);
}

} // namespace precice::time
