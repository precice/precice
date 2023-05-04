#include "time/Storage.hpp"
#include "math/differences.hpp"
#include "utils/assertion.hpp"

namespace precice::time {

const double Storage::WINDOW_START = 0.0;

const double Storage::WINDOW_END = 1.0;

Storage::Storage()
    : _sampleStorage{}, _extrapolationOrder{0}
{
}

void Storage::initialize(time::Sample sample)
{
  _sampleStorage.emplace_back(Stample{WINDOW_START, sample});
  _sampleStorage.emplace_back(Stample{WINDOW_END, sample});
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

void Storage::setExtrapolationOrder(int extrapolationOrder)
{
  _extrapolationOrder = extrapolationOrder;
}

double Storage::maxStoredNormalizedDt() const
{
  if (_sampleStorage.size() == 0) {
    return -1; // invalid return
  } else {
    return _sampleStorage.back().timestamp;
  }
}

int Storage::nTimes() const
{
  return _sampleStorage.size();
}

int Storage::nDofs() const
{
  PRECICE_ASSERT(_sampleStorage.size() > 0);
  return _sampleStorage[0].sample.values.size();
}

void Storage::move()
{
  PRECICE_ASSERT(nTimes() > 0);
  auto sampleAtBeginning = getSampleAtEnd();
  auto sampleAtEnd       = computeExtrapolation();
  _sampleStorage.clear();
  _sampleStorage.emplace_back(time::Stample{WINDOW_START, sampleAtBeginning});
  _sampleStorage.emplace_back(time::Stample{WINDOW_END, sampleAtEnd});
}

void Storage::clear()
{
  PRECICE_ASSERT(nTimes() > 0, "Storage does not contain any data!");
  Stample keep = _sampleStorage.front();
  PRECICE_ASSERT(math::equals(keep.timestamp, time::Storage::WINDOW_START));
  _sampleStorage.clear();
  _sampleStorage.emplace_back(keep);
}

void Storage::clearAll()
{
  _sampleStorage.clear();
}

Eigen::VectorXd Storage::getValuesAtOrAfter(double before) const
{
  auto sample = std::find_if(_sampleStorage.begin(), _sampleStorage.end(), [&before](const auto &s) { return math::greaterEquals(s.timestamp, before); });
  PRECICE_ASSERT(sample != _sampleStorage.end(), "no values found!");

  return sample->sample.values;
}

Eigen::VectorXd Storage::getTimes() const
{
  auto times = Eigen::VectorXd(nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i] = _sampleStorage[i].timestamp;
  }
  return times;
}

const std::vector<Stample> &Storage::getStamples() const
{
  PRECICE_DEBUG("Storage::getStamples()");
  return _sampleStorage;
}

std::pair<Eigen::VectorXd, Eigen::MatrixXd> Storage::getTimesAndValues() const
{
  auto times  = Eigen::VectorXd(nTimes());
  auto values = Eigen::MatrixXd(nDofs(), nTimes());
  for (int i = 0; i < times.size(); i++) {
    times[i]      = _sampleStorage[i].timestamp;
    values.col(i) = _sampleStorage[i].sample.values;
  }
  return std::make_pair(times, values);
}

time::Sample Storage::computeExtrapolation()
{
  if (_extrapolationOrder == 0) {
    return getSampleAtEnd(); // use values at end of window as initial guess for next
  } else if (_extrapolationOrder == 1) {
    auto s0 = getSampleAtBeginning();
    auto s1 = getSampleAtEnd();
    return time::Sample{2 * s1.values - s0.values, 2 * s1.gradient - s1.gradient}; // use linear extrapolation from window at beginning and end of window.
  }
  PRECICE_UNREACHABLE("Invalid _extrapolationOrder")
}

time::Sample Storage::getSampleAtBeginning()
{
  return _sampleStorage.front().sample;
}

time::Sample Storage::getSampleAtEnd()
{
  return _sampleStorage.back().sample;
}

} // namespace precice::time
