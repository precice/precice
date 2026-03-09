#pragma once

#include <chrono>
#include <vector>

namespace precice::impl {

class SolverImbalance {
  using Clock = std::chrono::steady_clock;

public:
  SolverImbalance() = default;
  double getSolverTimeToAdvance();
  void   startSolver();
  void   stopSolver(double solverTimeStepSize);
  void   computeSolverImbalance(const std::vector<double> &solverAdvanceTimes);
  void   reset();
  double getImbalance();
  double getImbalanceFactor();

private:
  std::vector<double>     _participant_dt;
  mutable logging::Logger _log{"impl::SolverImbalance"};
  std::vector<double>     _solver_advance_time;
  Clock::time_point       _start_time;
  double                  _solver_time_to_advance;
  double                  _imbalance;
  double                  _imbalanceFactor;

  enum class State : bool {
    STOPPED = false,
    RUNNING = true
  };

  State _state = State::STOPPED;
};

} // namespace precice::impl