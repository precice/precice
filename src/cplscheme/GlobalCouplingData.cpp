#include "cplscheme/GlobalCouplingData.hpp"

#include <utility>

#include "math/differences.hpp" // for math::equals(...)
#include "mesh/Data.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice::cplscheme {

GlobalCouplingData::GlobalCouplingData(
    mesh::PtrData data,
    bool          requiresInitialization,
    int           extrapolationOrder)
    : requiresInitialization(requiresInitialization),
      _data(std::move(data))
{
  PRECICE_ASSERT(_data != nullptr);
  /// Lazy allocation of _previousIteration.gradient: only used in case the corresponding data has gradients
  _previousIteration = time::Sample{Eigen::VectorXd::Zero(getSize())};
  timeStepsStorage().setExtrapolationOrder(extrapolationOrder);
}

int GlobalCouplingData::getDimensions() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->getDimensions();
}

int GlobalCouplingData::getSize() const
{
  return sample().values.size();
}

Eigen::VectorXd &GlobalCouplingData::values()
{
  return sample().values;
}

const Eigen::VectorXd &GlobalCouplingData::values() const
{
  return sample().values;
}

time::Storage &GlobalCouplingData::timeStepsStorage()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->timeStepsStorage();
}

const time::Storage &GlobalCouplingData::timeStepsStorage() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->timeStepsStorage();
}

void GlobalCouplingData::setSampleAtTime(double time, time::Sample sample)
{
  this->sample() = sample; // @todo at some point we should not need this anymore, when mapping, acceleration ... directly work on _timeStepsStorage
  timeStepsStorage().setSampleAtTime(time, sample);
}

void GlobalCouplingData::storeIteration()
{
  const auto stamples = this->stamples();
  PRECICE_ASSERT(stamples.size() > 0);
  this->sample()     = stamples.back().sample;
  _previousIteration = this->sample();
}

const Eigen::VectorXd GlobalCouplingData::previousIteration() const
{
  return _previousIteration.values;
}

int GlobalCouplingData::getPreviousIterationSize() const
{
  return _previousIteration.values.size();
}

int GlobalCouplingData::getDataID()
{
  return _data->getID();
}

std::string GlobalCouplingData::getDataName()
{
  return _data->getName();
}

void GlobalCouplingData::moveToNextWindow()
{
  if (this->timeStepsStorage().stamples().size() > 0) {
    this->timeStepsStorage().move();
    auto atEnd = this->timeStepsStorage().stamples().back();
    PRECICE_ASSERT(math::equals(atEnd.timestamp, time::Storage::WINDOW_END));
    _data->sample() = atEnd.sample;
  }
}

time::Sample &GlobalCouplingData::sample()
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->sample();
}

const time::Sample &GlobalCouplingData::sample() const
{
  PRECICE_ASSERT(_data != nullptr);
  return _data->sample();
}

} // namespace precice::cplscheme
