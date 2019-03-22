#pragma once

#include "../CouplingData.hpp"
#include "ConvergenceMeasure.hpp"
#include "logging/Logger.hpp"
#include "math/math.hpp"
#include "utils/MasterSlave.hpp"

namespace precice
{
namespace cplscheme
{
namespace tests
{
class RelativeConvergenceMeasureTest;
}
}
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice
{
namespace cplscheme
{
namespace impl
{

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
class RelativeConvergenceMeasure : public ConvergenceMeasure
{
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
      const Eigen::VectorXd &newValues,
      const Eigen::VectorXd &designSpecification)
  {
    /*
     std::cout<<"\n-------\n";
     std::cout<<"   old val: \n"<<oldValues<<'\n';
     std::cout<<"   new val: \n"<<newValues<<"\n\n";
     std::cout<<"   design spec: \n"<<designSpecification<<"\n\n";
     std::cout<<"-------\n\n";
*/

    _normDiff      = utils::MasterSlave::l2norm((newValues - oldValues) - designSpecification);
    _norm          = utils::MasterSlave::l2norm(newValues + designSpecification);
    _isConvergence = _normDiff <= _norm * _convergenceLimitPercent;
    //      INFO("Relative convergence measure: "
    //                    << "two-norm differences = " << normDiff
    //                    << ", convergence limit = "
    //                    << normNew * _convergenceLimitPercent
    //                    << ", convergence = " << _isConvergence );
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
    os << "two-norm diff = " << _normDiff;
    os << ", relative limit = " << _norm * _convergenceLimitPercent;
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

private:
  logging::Logger _log{"cplscheme::RelativeConvergenceMeasure"};

  double _convergenceLimitPercent;

  double _normDiff = 0;

  double _norm = 0;

  bool _isConvergence = false;
};
}
}
} // namespace precice, cplscheme, impl
