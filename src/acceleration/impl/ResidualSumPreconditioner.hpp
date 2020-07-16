#pragma once

#include <Eigen/Core>
#include <stddef.h>
#include <string>
#include <vector>
#include "acceleration/impl/Preconditioner.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace acceleration {
namespace impl {

/**
 * @brief Preconditioner that uses the residuals of all iterations of the current timestep summed up to scale the quasi-Newton system.
 * This is somewhat similar to what is done in the Marks and Luke paper.
 */
class ResidualSumPreconditioner : public Preconditioner {
public:
  ResidualSumPreconditioner(int maxNonConstTimesteps);
  /**
   * @brief Destructor, empty.
   */
  virtual ~ResidualSumPreconditioner() {}

  virtual void initialize(std::vector<size_t> &svs);

private:
  /**
   * @brief Update the scaling after every FSI iteration.
   *
   * @param[in] timestepComplete True if this FSI iteration also completed a timestep
   */
  virtual void _update_(bool timestepComplete, const Eigen::VectorXd &oldValues, const Eigen::VectorXd &res);

  logging::Logger _log{"acceleration::ResidualSumPreconditioner"};

  std::vector<double> _residualSum;
};

} // namespace impl
} // namespace acceleration
} // namespace precice
