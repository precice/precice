#pragma once

#include <Eigen/Core>
#include <ostream>
#include <string>
#include "ConvergenceMeasure.hpp"
#include "cplscheme/CouplingData.hpp"
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

class MinIterationConvergenceMeasure : public ConvergenceMeasure {
public:
  explicit MinIterationConvergenceMeasure(int minimumIterationCount);

  virtual ~MinIterationConvergenceMeasure() {}

  virtual void newMeasurementSeries();

  virtual void measure(
      const Eigen::VectorXd &oldValues,
      const Eigen::VectorXd &newValues)
  {
    PRECICE_TRACE();
    _currentIteration++;
    _isConvergence = (_minimumIterationCount <= _currentIteration);
    PRECICE_DEBUG("Iteration number = " << _currentIteration << ", convergence = " << _isConvergence);
  }

  virtual bool isConvergence() const
  {
    return _isConvergence;
  }

  virtual std::string printState()
  {
    std::ostringstream os;
    os << "min iteration convergence measure: ";
    os << "#it = " << _currentIteration << " of " << _minimumIterationCount;
    os << ", conv = ";
    if (_isConvergence)
      os << "true";
    else
      os << "false";
    return os.str();
  }

private:
  logging::Logger _log{"cplscheme::MinIterationConvergenceMeasure"};

  int _minimumIterationCount;

  int _currentIteration = 0;

  bool _isConvergence = false;
};
} // namespace impl
} // namespace cplscheme
} // namespace precice
