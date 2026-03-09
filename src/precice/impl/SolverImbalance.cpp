#include "precice/impl/SolverImbalance.hpp"
#include <algorithm>
#include <numeric>
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"

namespace precice::impl {

double SolverImbalance::getSolverTimeToAdvance()
{
  PRECICE_ASSERT(_participant_dt.size() == _solver_advance_time.size());
  double sum_pdt          = std::accumulate(_participant_dt.begin(), _participant_dt.end(), 0.0);
  double sum_sat          = std::accumulate(_solver_advance_time.begin(), _solver_advance_time.end(), 0.0);
  _solver_time_to_advance = sum_sat / sum_pdt;
  return _solver_time_to_advance;
}

void SolverImbalance::startSolver()
{
  PRECICE_ASSERT(_state == State::STOPPED);
  _state      = State::RUNNING;
  _start_time = Clock::now();
}

void SolverImbalance::stopSolver(double solverTimeStepSize)
{
  if (_state == State::STOPPED) { // first time step
    return;
  }
  _state                                = State::STOPPED;
  auto                     current_time = Clock::now();
  std::chrono::nanoseconds advance_time = current_time - _start_time;
  _solver_advance_time.push_back(advance_time.count() * 1e-9);
  _participant_dt.push_back(solverTimeStepSize);
  PRECICE_INFO("solver advance time: {}", _solver_advance_time.back());
  PRECICE_INFO("solver dt: {}", solverTimeStepSize);
  PRECICE_INFO("time per time step: {}", _solver_advance_time.back() / solverTimeStepSize);
}

void SolverImbalance::computeSolverImbalance(const std::vector<double> &solverAdvanceTimes)
{
  double mean_sat  = std::accumulate(solverAdvanceTimes.begin(), solverAdvanceTimes.end(), 0.0);
  mean_sat         = mean_sat / solverAdvanceTimes.size();
  auto max_time    = std::max_element(solverAdvanceTimes.begin(), solverAdvanceTimes.end());
  _imbalance       = *max_time / mean_sat;
  _imbalanceFactor = *max_time / _solver_time_to_advance;
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
  _participant_dt.clear();
  _solver_advance_time.clear();
}

} // namespace precice::impl