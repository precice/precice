#pragma once

#include <limits>
#include "../CouplingData.hpp"
#include "ConvergenceMeasure.hpp"
#include "logging/Logger.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

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

  virtual ~ResidualRelativeConvergenceMeasure(){};

  virtual void newMeasurementSeries()
  {
    _isConvergence     = false;
    _isFirstIteration  = true;
    _normFirstResidual = std::numeric_limits<double>::max();
  }

  virtual void measure(
      const Eigen::VectorXd &oldValues,
      const Eigen::VectorXd &newValues)
  {
    _normDiff = utils::MasterSlave::l2norm(newValues - oldValues);
    if (_isFirstIteration) {
      _normFirstResidual = _normDiff;
      _isFirstIteration  = false;
    }
    _isConvergence = _normDiff < _normFirstResidual * _convergenceLimitPercent;
    //      PRECICE_INFO("Residual Relative convergence measure: "
    //                    << "two-norm differences = " << normDiff
    //                    << ", convergence limit = "
    //                    << _normFirstResidual * _convergenceLimitPercent
    //                    << ", convergence = " << _isConvergence );
  }

  virtual bool isConvergence() const
  {
    return _isConvergence;
  }

  /// Adds current convergence information to output stream.
  virtual std::string printState()
  {
    std::ostringstream os;
    os << "residual relative convergence measure: ";
    os << "two-norm diff = " << _normDiff;
    os << ", relative limit = " << _normFirstResidual * _convergenceLimitPercent;
    os << ", conv = ";
    if (_isConvergence)
      os << "true";
    else
      os << "false";
    return os.str();
  }

private:
  logging::Logger _log{"cplscheme::ResidualRelativeConvergenceMeasure"};

  double _convergenceLimitPercent;

  bool _isFirstIteration = true;

  double _normFirstResidual = std::numeric_limits<double>::max();

  double _normDiff = 0;

  bool _isConvergence = false;
};
} // namespace impl
} // namespace cplscheme
} // namespace precice
