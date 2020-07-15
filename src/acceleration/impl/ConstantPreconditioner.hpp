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

/// Preconditioner that uses the constant user-defined factors to scale the quasi-Newton system.
class ConstantPreconditioner : public Preconditioner {
public:
  explicit ConstantPreconditioner(std::vector<double> factors);

  /**
   * @brief Destructor, empty.
   */
  virtual ~ConstantPreconditioner() {}

  virtual void initialize(std::vector<size_t> &svs);

private:
  /**
   * @brief Update the scaling after every FSI iteration.
   *
   * @param[in] timestepComplete True if this FSI iteration also completed a timestep
   */
  virtual void _update_(bool timestepComplete, const Eigen::VectorXd &oldValues, const Eigen::VectorXd &res);

  logging::Logger _log{"acceleration::ConstantPreconditioner"};

  /// Constant user-defined factors to scale the quasi-Newton system
  std::vector<double> _factors;
};

} // namespace impl
} // namespace acceleration
} // namespace precice
