#pragma once

#include <Eigen/Core>
#include <stddef.h>
#include <string>
#include <vector>
#include "acceleration/impl/Preconditioner.hpp"
#include "logging/Logger.hpp"

namespace precice::acceleration::impl {

/**
 * @brief Preconditioner that uses the residuals of all iterations of the current time window summed up to scale the quasi-Newton system.
 * This is somewhat similar to what is done in the Marks and Luke paper.
 */
class ResidualSumPreconditioner : public Preconditioner {
public:
  ResidualSumPreconditioner(int maxNonConstTimeWindows);
  /**
   * @brief Destructor, empty.
   */
  ~ResidualSumPreconditioner() override = default;

  void initialize(std::vector<size_t> &svs) override;

private:
  /**
   * @brief Update the scaling after every FSI iteration.
   *
   * @param[in] timeWindowComplete True if this FSI iteration also completed a time window
   */
  void _update_(bool timeWindowComplete, const Eigen::VectorXd &oldValues, const Eigen::VectorXd &res) override;

  logging::Logger _log{"acceleration::ResidualSumPreconditioner"};

  std::vector<double> _residualSum;
};

} // namespace precice::acceleration::impl
