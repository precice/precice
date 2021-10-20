#pragma once

#include <Eigen/Core>
#include <string>
#include "acceleration/impl/Preconditioner.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace acceleration {
namespace impl {

/**
 * @brief Preconditioner that uses the recent residual to scale the quasi-Newton system.
 */
class ResidualPreconditioner : public Preconditioner {
public:
  ResidualPreconditioner(
      int maxNonConstTimeWindows);

  /**
   * @brief Destructor, empty.
   */
  virtual ~ResidualPreconditioner() {}

private:
  /**
    * @brief Update the scaling after every FSI iteration.
    *
    * @param[in] timeWindowComplete True if this FSI iteration also completed a time window
    */
  virtual void _update_(bool                   timeWindowComplete,
                        const Eigen::VectorXd &oldValues,
                        const Eigen::VectorXd &res);

  logging::Logger _log{"acceleration::ResidualPreconditioner"};
};

} // namespace impl
} // namespace acceleration
} // namespace precice
