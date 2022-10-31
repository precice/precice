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

Eigen::VectorXd Storage::getValueAtTime(double time)
{
  for (auto &sample : _sampleStorage) {
    if (math::equals(sample.first, time)) {
      return sample.second;
    }
  }
  PRECICE_ASSERT(false, "no value found!", time);
}

void Storage::setValueAtTime(double time, Eigen::VectorXd value)
{
  PRECICE_ASSERT(math::greater(time, WINDOW_START), "Setting value outside of valid range!");
  PRECICE_ASSERT(math::smallerEquals(time, WINDOW_END), "Sampling outside of valid range!");
  PRECICE_ASSERT(math::smaller(maxStoredNormalizedDt(), time), maxStoredNormalizedDt(), time, "Trying to overwrite existing values or to write values with a time that is too small. Please use clear(), if you want to reset the storage.");
  _sampleStorage.emplace_back(std::make_pair(time, value));
}

void Storage::overrideDataAtEndWindowTime(Eigen::VectorXd data)
{
  if (_sampleStorage.size() == 0) {
    _sampleStorage.emplace_back(std::make_pair(WINDOW_END, data));
  } else {
    PRECICE_ASSERT(math::equals(_sampleStorage.back().first, WINDOW_END), "Unexpected!", _sampleStorage.back().first);
    _sampleStorage.back().second = data;
  }
}

double Storage::maxStoredNormalizedDt()
{
  if (_sampleStorage.size() == 0) {
    return 0; // @todo better return something that is clearly invalid or raise an error.
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
  auto initialGuess = _sampleStorage.back().second; // use value at end of window as initial guess for next
  _sampleStorage.clear();
  initialize(initialGuess);
}

void Storage::clear(bool keepZero)
{
  Eigen::VectorXd keep;
  if (keepZero) {
    keep = _sampleStorage.front().second; // we keep data at _storageDict[0.0]
  }
  _sampleStorage.clear();
  if (keepZero) {
    _sampleStorage.emplace_back(std::make_pair(WINDOW_START, keep));
  }
}

Eigen::VectorXd Storage::getValueAtOrAfter(double before)
{
  for (auto &sample : _sampleStorage) {
    if (math::greaterEquals(sample.first, before)) {
      return sample.second;
    }
  }
  PRECICE_ASSERT(false, "no value found!", before);
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

} // namespace precice::time
