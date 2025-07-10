#pragma once

#include <Eigen/Core>
#include <iomanip>
#include <limits>
#include <ostream>
#include <string>
#include "../CouplingData.hpp"
#include "ConvergenceMeasure.hpp"
#include "logging/Logger.hpp"
#include "math/differences.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme::impl {

/**
 * @brief Measures the convergence from an old data set to a new one.
 *
 * The convergence is evaluated by looking at the two norm of the differences
 * between each data value from the new and old data set. If the two norm is
 * equal or below a given percentage of the norm of the old data set,
 * convergence is achieved.
 *
 * For a description of how to perform the measurement, see class
 * ConvergenceMeasure.
 */
class ResidualRelativeConvergenceMeasure : public ConvergenceMeasure {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] convergenceLimitPercent
   *        Limit to define convergence relative to the norm of the current
   *        new dataset. Has to be in $] 0 ; 1 ]$.
   */
  explicit ResidualRelativeConvergenceMeasure(double convergenceLimitPercent);

  ~ResidualRelativeConvergenceMeasure() override = default;

  void newMeasurementSeries() override
  {
    _isConvergence     = false;
    _isFirstIteration  = true;
    _normFirstResidual = std::numeric_limits<double>::max();
  }

  void measure(
      const Eigen::VectorXd &oldValues,
      const Eigen::VectorXd &newValues) override
  {
    PRECICE_ASSERT(oldValues.size() == newValues.size());
    _normDiff = utils::IntraComm::l2norm(newValues - oldValues);
    if (_isFirstIteration) {
      _normFirstResidual = _normDiff;
      _isFirstIteration  = false;
    }
    _isConvergence = _normDiff < _normFirstResidual * _convergenceLimitPercent;
  }

  bool isConvergence() const override
  {
    return _isConvergence;
  }

  /// Adds current convergence information to output stream.
  std::string printState(const std::string &dataName) override
  {
    std::ostringstream os;
    os << "residual relative convergence measure: ";
    os << "relative two-norm diff of data \"" << dataName << "\" = ";
    os << std::scientific << std::setprecision(2) << getNormResidual();
    os << ", limit = " << _convergenceLimitPercent;
    os << ", normalization = " << _normFirstResidual;
    os << ", conv = ";
    if (_isConvergence)
      os << "true";
    else
      os << "false";
    return os.str();
  }

  double getNormResidual() override
  {
    if (math::equals(_normFirstResidual, 0.))
      return std::numeric_limits<double>::infinity();
    else
      return _normDiff / _normFirstResidual;
  }

  std::string getAbbreviation() const override
  {
    return "Drop";
  }

private:
  logging::Logger _log{"cplscheme::ResidualRelativeConvergenceMeasure"};

  double _convergenceLimitPercent;

  bool _isFirstIteration = true;

  double _normFirstResidual = std::numeric_limits<double>::max();

  double _normDiff = 0;

  bool _isConvergence = false;
};
} // namespace precice::cplscheme::impl
