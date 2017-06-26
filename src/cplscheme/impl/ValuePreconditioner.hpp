#pragma once

#include "Preconditioner.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Preconditioner that uses the values from the previous timestep to scale the quasi-Newton system.
 */
class ValuePreconditioner : public Preconditioner
{
public:

  ValuePreconditioner(
      int maxNonConstTimesteps);
  /**
   * @brief Destructor, empty.
   */
  virtual ~ValuePreconditioner() {}


private:

  /**
   * @brief Update the scaling after every FSI iteration.
   *
   * @param timestepComplete [IN] True if this FSI iteration also completed a timestep
   */
  virtual void _update_(bool timestepComplete, const Eigen::VectorXd& oldValues, const Eigen::VectorXd& res);

  static logging::Logger _log;

  bool _firstTimestep;
};

}}} // namespace precice, cplscheme, impl
