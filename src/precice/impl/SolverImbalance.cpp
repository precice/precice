#include "precice/impl/SolverImbalance.hpp"
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
  PRECICE_INFO("starting solver imbalance");
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

} // namespace precice::impl