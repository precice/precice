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

void Storage::setValueAtTime(double time, Eigen::VectorXd value)
{
  PRECICE_ASSERT(math::greater(time, WINDOW_START), "Setting value outside of valid range!");
  PRECICE_ASSERT(math::smallerEquals(time, WINDOW_END), "Sampling outside of valid range!");
  PRECICE_ASSERT(math::smaller(maxStoredNormalizedDt(), time), "Trying to overwrite existing values or to write values with a time that is too small. Please use clear(), if you want to reset the storage.");
  _sampleStorage.emplace_back(std::make_pair(time, value));
}

double Storage::maxStoredNormalizedDt()
{
  PRECICE_ASSERT(_sampleStorage.size() > 0);
  return _sampleStorage.back().first;
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
  auto sample = std::find_if(_sampleStorage.begin(), _sampleStorage.end(), [&before](const auto &s) { return math::greaterEquals(s.first, before); });
  PRECICE_ASSERT(sample != _sampleStorage.end(), "no value found!");

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

} // namespace precice::time
