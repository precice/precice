#include "precice/impl/SolverImbalance.hpp"
#include <algorithm>
#include <numeric>
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"

namespace precice::impl {

double SolverImbalance::getSolverTimeToAdvance()
{
  _solverTimeToAdvance = _solverAdvanceTime / _simulatedTime;
  return _solverTimeToAdvance;
}

void SolverImbalance::startSolver()
{
  PRECICE_ASSERT(_state == State::STOPPED);
  _state     = State::RUNNING;
  _startTime = Clock::now();
}

void SolverImbalance::stopSolver(double solverTimeStepSize)
{
  if (_state == State::STOPPED) { // first time step
    return;
  }
  _state             = State::STOPPED;
  auto   currentTime = Clock::now();
  double advanceTime = (currentTime - _startTime).count() * 1e-9;
  _solverAdvanceTime += advanceTime;
  _simulatedTime += solverTimeStepSize;
}

void SolverImbalance::computeSolverImbalance(const std::vector<double> &timesToAdvance)
{
  double meanTimeToAdvance = std::accumulate(timesToAdvance.begin(), timesToAdvance.end(), 0.0);
  meanTimeToAdvance        = meanTimeToAdvance / timesToAdvance.size();
  auto maxTimeToAdvance    = std::max_element(timesToAdvance.begin(), timesToAdvance.end());
  _imbalance               = *maxTimeToAdvance / meanTimeToAdvance;
  _imbalanceFactor         = *maxTimeToAdvance / _solverTimeToAdvance;
}

double SolverImbalance::getImbalance()
{
  return _imbalance;
}

double SolverImbalance::getImbalanceFactor()
{
  return _imbalanceFactor;
}

void SolverImbalance::reset()
{
  PRECICE_DEBUG("resetting SolverImbalance");
  _simulatedTime     = 0.0;
  _solverAdvanceTime = 0.0;
}

} // namespace precice::impl
