#include "time/Storage.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

Storage::Storage()
    : _sampleStorage{}
{
}

void Storage::initialize(Eigen::VectorXd values)
{
  _sampleStorage.emplace_back(std::make_pair(0.0, values));
  _sampleStorage.emplace_back(std::make_pair(1.0, values));
}

int findTime(std::vector<std::pair<double, Eigen::VectorXd>> storage, double time)
{
  for (int i; i < storage.size(); ++i) {
    if (math::equals(storage[i].first, time)) {
      return i;
    }
  }
  return -1;
}

void Storage::setValueAtTime(double time, Eigen::VectorXd value)
{
  PRECICE_ASSERT(math::greater(time, 0.0), "Setting value outside of valid range!");
  PRECICE_ASSERT(math::smallerEquals(time, 1.0), "Sampling outside of valid range!");
  PRECICE_ASSERT(math::smaller(maxStoredNormalizedDt(), time));
  _sampleStorage.emplace_back(std::make_pair(time, value));
}

double Storage::maxStoredNormalizedDt()
{
  return _sampleStorage.back().first;
}

int Storage::nTimes()
{
  return _sampleStorage.size();
}

int Storage::nDofs()
{
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
    _sampleStorage.emplace_back(std::make_pair(0.0, keep));
  }
}

Eigen::VectorXd Storage::getValueAtTimeAfter(double before)
{
  for (auto &sample : _sampleStorage) {
    if (math::greaterEquals(_sampleStorage.first, before)) {
      return _sampleStorage.second;
    }
  }
  PRECICE_ASSERT(false, "no value found!");
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
