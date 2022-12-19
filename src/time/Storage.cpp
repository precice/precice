#include "time/Storage.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

const double Storage::WINDOW_START = 0.0;

const double Storage::WINDOW_END = 1.0;

Storage::Storage(int extrapolationOrder)
    : _sampleStorage{}, _extrapolationOrder(extrapolationOrder)
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

double Storage::maxStoredNormalizedDt()
{
  if (_sampleStorage.size() == 0) {
    return -1; // invalid return
  } else {
    return _sampleStorage.back().first;
  }
}

int Storage::nTimes()
{
  return _sampleStorage.size();
}

int Storage::nDofs()
{
  PRECICE_ASSERT(_sampleStorage.size() > 0);
  return _sampleStorage[0].second.size();
}

void Storage::move()
{
  PRECICE_ASSERT(nTimes() > 0);
  auto valuesAtBeginning = getValuesAtEnd();
  auto valuesAtEnd       = computeExtrapolation();
  _sampleStorage.clear();
  _sampleStorage.emplace_back(std::make_pair(WINDOW_START, valuesAtBeginning));
  _sampleStorage.emplace_back(std::make_pair(WINDOW_END, valuesAtEnd));
}

void Storage::clear()
{
  PRECICE_ASSERT(nTimes() > 0, "Storage does not contain any data!");
  Eigen::VectorXd keep = getValuesAtBeginning(); // we keep data at _storageDict[0.0]
  _sampleStorage.clear();
  _sampleStorage.emplace_back(std::make_pair(WINDOW_START, keep));
}

Eigen::VectorXd Storage::getValuesAtOrAfter(double before)
{
  auto sample = std::find_if(_sampleStorage.begin(), _sampleStorage.end(), [&before](const auto &s) { return math::greaterEquals(s.first, before); });
  PRECICE_ASSERT(sample != _sampleStorage.end(), "no values found!");

  return sample->second;
}

Eigen::VectorXd Storage::getTimes()
{
  auto times = Eigen::VectorXd(nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i] = _sampleStorage[i].first;
  }
  return times;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> Storage::getTimesAndValues()
{
  auto times  = Eigen::VectorXd(nTimes());
  auto values = Eigen::MatrixXd(nDofs(), nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i]      = _sampleStorage[i].first;
    values.col(i) = _sampleStorage[i].second;
  }
  return std::make_pair(times, values);
}

Eigen::VectorXd Storage::computeExtrapolation()
{
  if (_extrapolationOrder == 0) {
    return getValuesAtEnd(); // use values at end of window as initial guess for next
  } else if (_extrapolationOrder == 1) {
    return 2 * getValuesAtEnd() - getValuesAtBeginning(); // use linear extrapolation from window at beginning and end of window.
  }
  PRECICE_UNREACHABLE("Invalid _extrapolationOrder")
}

Eigen::VectorXd Storage::getValuesAtBeginning()
{
  PRECICE_ASSERT(_sampleStorage.front().first == WINDOW_START, _sampleStorage.front().first);
  return _sampleStorage.front().second;
}

Eigen::VectorXd Storage::getValuesAtEnd()
{
  PRECICE_ASSERT(_sampleStorage.back().first == WINDOW_END, _sampleStorage.back().first);
  return _sampleStorage.back().second;
}

} // namespace precice::time
