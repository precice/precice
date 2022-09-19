#include "time/Storage.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

Storage::Storage()
    : _storageDict{} {}

void Storage::initialize(Eigen::VectorXd values)
{
  _storageDict[1.0] = Eigen::VectorXd(values);
  _storageDict[0.0] = Eigen::VectorXd(values);
}

Eigen::VectorXd Storage::getValueAtTime(double time)
{
  PRECICE_ASSERT(time >= 0, "Sampling outside of valid range!");
  PRECICE_ASSERT(time <= 1, "Sampling outside of valid range!");
  PRECICE_ASSERT(_storageDict.size() > 0);
  PRECICE_ASSERT(_storageDict.count(time) > 0, time);
  return _storageDict[time];
}

void Storage::setValueAtTime(double time, Eigen::VectorXd value)
{
  PRECICE_ASSERT(time > 0, "Setting value outside of valid range!");
  PRECICE_ASSERT(time <= 1, "Sampling outside of valid range!");
  _storageDict[time] = value;
}

double Storage::maxStoredNormalizedDt()
{
  double theMaxKey = -1 * std::numeric_limits<double>::infinity();
  for (auto pair : _storageDict) {
    if (pair.first > theMaxKey) {
      theMaxKey = pair.first;
    }
  }
  PRECICE_ASSERT(theMaxKey > -1 * std::numeric_limits<double>::infinity());
  return theMaxKey;
}

int Storage::nTimes()
{
  return _storageDict.size();
}

int Storage::nDofs()
{
  return _storageDict[0.0].size();
}

void Storage::move()
{
  PRECICE_ASSERT(nTimes() > 0);
  auto initialGuess = _storageDict[maxStoredNormalizedDt()]; // use value at end of window as initial guess for next
  _storageDict.clear();
  _storageDict[0.0] = Eigen::VectorXd(initialGuess);
  _storageDict[1.0] = Eigen::VectorXd(initialGuess); // initial guess is always constant extrapolation
}

void Storage::clear(bool keepZero)
{
  Eigen::VectorXd keep;
  if (keepZero) {
    keep = _storageDict[0.0]; // we keep data at _storageDict[0.0]
  }
  _storageDict.clear();
  if (keepZero) {
    _storageDict[0.0] = keep;
  }
}

double Storage::getClosestTimeAfter(double before)
{
  double directlyAfter = std::numeric_limits<double>::infinity();

  for (auto pair : _storageDict) {
    if (before <= pair.first && pair.first <= directlyAfter) { // pair.first is after "before" and earlier than current "directlyAfter"
      directlyAfter = pair.first;
    }
  }

  return directlyAfter;
}

Eigen::VectorXd Storage::getTimes()
{
  // create std::vector with all keys
  std::vector<double> keys;
  for (auto timeStep : _storageDict) {
    keys.push_back(timeStep.first);
  }

  // sort vector
  std::sort(keys.begin(), keys.end());

  // copy data into Eigen::VectorXd to return
  auto times = Eigen::VectorXd(keys.size());
  for (int i = 0; i < keys.size(); i++) {
    times[i] = keys[i];
  }
  return times;
}

} // namespace precice::time
