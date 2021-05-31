#pragma once

#include <Eigen/Core>
#include <limits>
#include <math.h>
#include <ostream>
#include <string>
#include "../CouplingData.hpp"
#include "ConvergenceMeasure.hpp"
#include "logging/Logger.hpp"
#include "math/differences.hpp"
#include "math/math.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {
namespace tests {
class RelativeConvergenceMeasureTest;
}
} // namespace cplscheme
} // namespace precice

// ----------------------------------------------------------- CLASS DEFINITION

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
class RelativeConvergenceMeasure : public ConvergenceMeasure {
public:
  /**
    * @brief Constructor.
    *
    * @param[in] convergenceLimitPercent
    *        Limit to define convergence relative to the norm of the current
    *        new dataset. Has to be in $] 0 ; 1 ]$.
    */
  RelativeConvergenceMeasure(double convergenceLimitPercent);

  virtual ~RelativeConvergenceMeasure(){};

  virtual void newMeasurementSeries()
  {
    _isConvergence = false;
  }

  virtual void measure(
      const Eigen::VectorXd &oldValues,
      const Eigen::VectorXd &newValues)
  {
    /*
     std::cout<<"\n-------\n";
     std::cout<<"   old val: \n"<<oldValues<<'\n';
     std::cout<<"   new val: \n"<<newValues<<"\n\n";
     std::cout<<"-------\n\n";
*/

    _normDiff      = utils::MasterSlave::l2norm(newValues - oldValues);
    _norm          = utils::MasterSlave::l2norm(newValues);
    _isConvergence = _normDiff <= _norm * _convergenceLimitPercent;
  }

  virtual bool isConvergence() const
  {
    return _isConvergence;
  }

  /**
    * @brief Adds current convergence information to output stream.
    */
  virtual std::string printState()
  {
    std::ostringstream os;
    os << "relative convergence measure: ";
    os << "relative two-norm diff = " << getNormResidual();
    os << ", limit = " << _convergenceLimitPercent;
    os << ", normalization = " << _norm;
    os << ", conv = ";
    if (_isConvergence)
      os << "true";
    else
      os << "false";
    return os.str();
  }

  virtual double getNormResidual()
  {
    if (math::equals(_norm, 0.))
      return std::numeric_limits<double>::infinity();
    else
      return _normDiff / _norm;
  }

  virtual std::string getAbbreviation() const
  {
    return "Rel";
  }

private:
  logging::Logger _log{"cplscheme::RelativeConvergenceMeasure"};

  double _convergenceLimitPercent;

  double _normDiff = 0;

  double _norm = 0;

  bool _isConvergence = false;
};
} // namespace impl
} // namespace cplscheme
} // namespace precice
