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
  void   computeSolverImbalance(const std::vector<double> &timesToAdvance);
  void   reset();
  double getImbalance();
  double getImbalanceFactor();

private:
  mutable logging::Logger _log{"impl::SolverImbalance"};
  double                  _solverAdvanceTime = 0.0;
  double                  _simulatedTime     = 0.0;
  Clock::time_point       _startTime;
  double                  _solverTimeToAdvance;
  double                  _imbalance;
  double                  _imbalanceFactor;

  enum class State : bool {
    STOPPED = false,
    RUNNING = true
  };

  State _state = State::STOPPED;
};

} // namespace precice::impl
