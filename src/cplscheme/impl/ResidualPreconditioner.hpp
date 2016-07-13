#pragma once

#include "Preconditioner.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "utils/Globals.hpp"
#include "Eigen/Dense"
#include "../SharedPointer.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Preconditioner that uses the recent residual to scale the quasi-Newton system.
 */
class ResidualPreconditioner : public Preconditioner
{
public:

  ResidualPreconditioner(
      std::vector<int> dimensions,
      int maxNonConstTimesteps);

  /**
   * @brief Destructor, empty.
   */
  virtual ~ResidualPreconditioner() {}


private:

  /**
    * @brief Update the scaling after every FSI iteration.
    *
    * @param timestepComplete [IN] True if this FSI iteration also completed a timestep
    */
   virtual void _update_(bool timestepComplete, const Eigen::VectorXd& oldValues, const Eigen::VectorXd& res);


  static tarch::logging::Log _log;
};

}}} // namespace precice, cplscheme, impl
