#pragma once

#include <Eigen/Core>
#include <iomanip>
#include <ostream>
#include <string>
#include "ConvergenceMeasure.hpp"
#include "logging/Logger.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme::tests {
class AbsoluteConvergenceMeasureTest;
}

namespace precice::cplscheme::impl {

/**
 * @brief Measures the convergence from an old data set to a new one.
 *
 * The convergence is evaluated by looking at the two norm of the differences
 * between each data value from the new and old data set. If the two norm is
 * equal or below a given limit, convergence is achieved.
 *
 * For a description of how to perform the measurement, see class
 * ConvergenceMeasure.
 */
class AbsoluteConvergenceMeasure : public ConvergenceMeasure {
public:
  explicit AbsoluteConvergenceMeasure(double convergenceLimit);

  ~AbsoluteConvergenceMeasure() override = default;

  void newMeasurementSeries() override
  {
    _isConvergence = false;
  }

  void measure(
      const Eigen::VectorXd &oldValues,
      const Eigen::VectorXd &newValues) override
  {
    PRECICE_ASSERT(oldValues.size() == newValues.size());
    _normDiff      = utils::IntraComm::l2norm(newValues - oldValues);
    _isConvergence = _normDiff <= _convergenceLimit;
  }

  bool isConvergence() const override
  {
    return _isConvergence;
  }

  /// Adds current convergence information to output stream.
  std::string printState(const std::string &dataName) override
  {
    std::ostringstream os;
    os << "absolute convergence measure: ";
    os << "two-norm diff of data \"" << dataName << "\" = ";
    os << std::scientific << std::setprecision(2) << _normDiff;
    os << ", limit = " << _convergenceLimit;
    os << ", conv = ";
    if (_isConvergence)
      os << "true";
    else
      os << "false";
    return os.str();
  }

  double getNormResidual() override
  {
    return _normDiff;
  }

  std::string getAbbreviation() const override
  {
    return "Abs";
  }

private:
  logging::Logger _log{"cplscheme::AbsoluteConvergenceMeasure"};

  double _convergenceLimit;

  double _normDiff = 0;

  bool _isConvergence = false;
};
} // namespace precice::cplscheme::impl
