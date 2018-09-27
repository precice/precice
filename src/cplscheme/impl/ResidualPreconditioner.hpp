#pragma once

#include <Eigen/Core>
#include "Preconditioner.hpp"

namespace precice
{
namespace cplscheme
{
namespace impl
{

/**
 * @brief Preconditioner that uses the recent residual to scale the quasi-Newton system.
 */
class ResidualPreconditioner : public Preconditioner
{
public:
  ResidualPreconditioner(
      int maxNonConstTimesteps);

  /**
   * @brief Destructor, empty.
   */
  virtual ~ResidualPreconditioner() {}

private:
  /**
    * @brief Update the scaling after every FSI iteration.
    *
    * @param[in] timestepComplete True if this FSI iteration also completed a timestep
    */
  virtual void _update_(bool timestepComplete,
                        const Eigen::VectorXd &oldValues,
                        const Eigen::VectorXd &res);

  logging::Logger _log{"cplscheme::ResidualPreconditioner"};
};
}
}
} // namespace precice, cplscheme, impl
