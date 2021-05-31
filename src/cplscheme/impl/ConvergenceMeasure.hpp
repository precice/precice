#pragma once

#include <Eigen/Core>

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Interface for measures checking the convergence of a series of datasets.
 *
 * A measurement involves two states of the data set: an old state and a new
 * state. Typically, the states corresponds to timestep $t$ and $t+1$. The
 * subclasses of ConvergenceMeasure define how exactly convergence is measured.
 *
 * A measure has to be used in the following way:
 * -# create the measure object (a subclass of ConvergenceMeasure)
 * -# call newMeasurementSeries() for one set of iterations
 * -# call measure() for convergence measurement
 * -# retrieve the convergence status via isConvergence()
 */
class ConvergenceMeasure {
public:
  /// Destructor, empty.
  virtual ~ConvergenceMeasure() {}

  /// To be called when a new meas. series (iteration process) starts.
  virtual void newMeasurementSeries() = 0;

  /**
   * @brief Performs convergence measurement.
   *
   * @param[in] oldValues Old iterate values.
   * @param[in] newValues New iterate values.
   */
  virtual void measure(
      const Eigen::VectorXd &oldValues,
      const Eigen::VectorXd &newValues) = 0;

  /// Returns true, if the last measurement indicates convergence.
  virtual bool isConvergence() const = 0;

  /// Adds current convergence information to output stream.
  virtual std::string printState() = 0;

  /// Returns the l2-norm of the coupling residuum
  virtual double getNormResidual()
  {
    return 0;
  }

  /// Returns an abreviation of the name of the measure for the log file headers
  virtual std::string getAbbreviation() const
  {
    return "";
  }
};
} // namespace impl
} // namespace cplscheme
} // namespace precice
