#pragma once

#include <Eigen/Core>
#include <string>
#include "acceleration/impl/Preconditioner.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace acceleration {
namespace impl {

/// Preconditioner that uses the values from the previous timestep to scale the quasi-Newton system.
class ValuePreconditioner : public Preconditioner {
public:
  ValuePreconditioner(
      int maxNonConstTimesteps);
  /**
   * @brief Destructor, empty.
   */
  virtual ~ValuePreconditioner() {}

private:
  logging::Logger _log{"acceleration::ValuePreconditioner"};

  /**
   * @brief Update the scaling after every FSI iteration.
   *
   * @param[in] timestepComplete True if this FSI iteration also completed a timestep
   */
  virtual void _update_(bool timestepComplete, const Eigen::VectorXd &oldValues, const Eigen::VectorXd &res);

  bool _firstTimestep = true;
};

} // namespace impl
} // namespace acceleration
} // namespace precice
