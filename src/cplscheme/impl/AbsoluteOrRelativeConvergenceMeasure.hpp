#pragma once

#include <Eigen/Core>
#include <iomanip>
#include <limits>
#include <math.h>
#include <ostream>
#include <string>
#include "../CouplingData.hpp"
#include "ConvergenceMeasure.hpp"
#include "logging/Logger.hpp"
#include "math/differences.hpp"
#include "math/math.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::cplscheme::tests {
class AbsoluteOrRelativeConvergenceMeasureTest;
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice::cplscheme::impl {

/**
 * @brief Measures the convergence from an old data set to a new one.
 *
 * The convergence is evaluated by looking at the two norm of the differences
 * between each data value from the new and old data set. If the two norm is
 * equal or below a given absolute value or a given percentage of the norm of
 * the new data set,
 * convergence is achieved.
 *
 * For a description of how to perform the measurement, see class
 * ConvergenceMeasure.
 */
class AbsoluteOrRelativeConvergenceMeasure : public ConvergenceMeasure {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] absLimit
   *        Limit to define absolute convergence to the norm of the current
   *        new dataset. Has to be larger than zero.
   * @param[in] relLimit
   *        Limit to define convergence relative to the norm of the current
   *        new dataset. Has to be in $] 0 ; 1 ]$.
   */
  AbsoluteOrRelativeConvergenceMeasure(double absLimit, double relLimit);

  ~AbsoluteOrRelativeConvergenceMeasure() override = default;

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
    _norm          = utils::IntraComm::l2norm(newValues);
    _isConvergence = (_normDiff <= _norm * _convergenceLimitPercent) or (_normDiff <= _convergenceLimit);
  }

  bool isConvergence() const override
  {
    return _isConvergence;
  }

  /**
   * @brief Adds current convergence information to output stream.
   */
  std::string printState(const std::string &dataName) override
  {
    std::ostringstream os;
    os << "absolute convergence measure: ";
    os << "two-norm diff of data \"" << dataName << "\" = ";
    os << std::scientific << std::setprecision(2) << getNormAbsResidual();
    os << ", limit = " << _convergenceLimit;

    os << ", relative convergence measure: ";
    os << "relative two-norm diff of data \"" << dataName << "\" = ";
    os << std::scientific << std::setprecision(2) << getNormRelResidual();
    os << ", limit = " << _convergenceLimitPercent;
    os << ", normalization = " << _norm;
    os << ", conv = ";
    if (_isConvergence)
      os << "true";
    else
      os << "false";
    return os.str();
  }

  virtual double getNormAbsResidual()
  {
    return _normDiff;
  }

  virtual double getNormRelResidual()
  {
    if (math::equals(_norm, 0.))
      return std::numeric_limits<double>::infinity();
    else
      return _normDiff / _norm;
  }

  std::string getAbbreviation() const override
  {
    return "AbsOrRel";
  }

private:
  logging::Logger _log{"cplscheme::AbsoluteOrRelativeConvergenceMeasure"};

  double _convergenceLimit;

  double _convergenceLimitPercent;

  double _normDiff = 0;

  double _norm = 0;

  bool _isConvergence = false;
};
} // namespace precice::cplscheme::impl
