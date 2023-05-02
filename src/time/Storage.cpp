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
  _sampleStorage.emplace_back(std::make_pair(WINDOW_START, values));
  _sampleStorage.emplace_back(std::make_pair(WINDOW_END, values));
}

Eigen::VectorXd Storage::getValuesAtTime(double time)
{
  for (auto &sample : _sampleStorage) {
    if (math::equals(sample.first, time)) {
      return sample.second;
    }
  }
  PRECICE_ASSERT(false, "no values found!", time);
}

void Storage::setValuesAtTime(double time, Eigen::VectorXd values, bool mustOverwriteExisting)
{
  PRECICE_ASSERT(math::smallerEquals(WINDOW_START, time), "Setting values outside of valid range!");
  PRECICE_ASSERT(math::smallerEquals(time, WINDOW_END), "Sampling outside of valid range!");
  if (!mustOverwriteExisting) {
    PRECICE_ASSERT(math::smaller(maxStoredNormalizedDt(), time), maxStoredNormalizedDt(), time, "Trying to overwrite existing values or to write values with a time that is too small. Please use clear(), if you want to reset the storage.");
    _sampleStorage.emplace_back(std::make_pair(time, values));
  } else {
    // check that key "time" exists.
    auto sample = std::find_if(_sampleStorage.begin(), _sampleStorage.end(), [&time](const auto &s) { return math::equals(s.first, time); });
    PRECICE_ASSERT(sample != _sampleStorage.end(), time, "Key does not exist, cannot overwrite values.");
    // overwrite values at "time"
    for (auto &sample : _sampleStorage) {
      if (math::equals(sample.first, time)) {
        sample.second = values;
        return;
      }
    }
    PRECICE_ASSERT(false, "unreachable!");
  }
}

double Storage::maxStoredNormalizedDt() const
{
  if (_sampleStorage.size() == 0) {
    return -1; // invalid return
  } else {
    return _sampleStorage.back().first;
  }
}

int Storage::nTimes() const
{
  return _sampleStorage.size();
}

int Storage::nDofs() const
{
  PRECICE_ASSERT(_sampleStorage.size() > 0);
  return _sampleStorage[0].second.size();
}

void Storage::move()
{
  PRECICE_ASSERT(nTimes() > 0);
  auto initialGuess = _sampleStorage.back().second; // use values at end of window as initial guess for next
  _sampleStorage.clear();
  initialize(initialGuess);
}

void Storage::clear(bool keepWindowStart)
{
  if (keepWindowStart) {
    PRECICE_ASSERT(nTimes() > 0, "Storage does not contain any data!");
    Eigen::VectorXd keep = _sampleStorage.front().second; // we keep data at _storageDict[0.0]
    _sampleStorage.clear();
    _sampleStorage.emplace_back(std::make_pair(WINDOW_START, keep));
  } else {
    _sampleStorage.clear();
  }
}

Eigen::VectorXd Storage::getValuesAtOrAfter(double before) const
{
  auto sample = std::find_if(_sampleStorage.begin(), _sampleStorage.end(), [&before](const auto &s) { return math::greaterEquals(s.first, before); });
  PRECICE_ASSERT(sample != _sampleStorage.end(), "no values found!");

  return sample->second;
}

Eigen::VectorXd Storage::getTimes() const
{
  auto times = Eigen::VectorXd(nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i] = _sampleStorage[i].first;
  }
  return times;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> Storage::getTimesAndValues() const
{
  auto times  = Eigen::VectorXd(nTimes());
  auto values = Eigen::MatrixXd(nDofs(), nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i]      = _sampleStorage[i].first;
    values.col(i) = _sampleStorage[i].second;
  }
  return std::make_pair(times, values);
}

} // namespace precice::time
