// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#pragma once

#include "Preconditioner.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "utils/Globals.hpp"
#include "../SharedPointer.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Preconditioner that uses the constant user-defined factors to scale the quasi-Newton system.
 */
class ConstantPreconditioner : public Preconditioner
{
public:

  ConstantPreconditioner(
      std::vector<int> dimensions,
      std::vector<double> factors
  );
  /**
   * @brief Destructor, empty.
   */
  virtual ~ConstantPreconditioner() {}

  virtual void initialize(int N);

private:

  /**
   * @brief Update the scaling after every FSI iteration.
   *
   * @param timestepComplete [IN] True if this FSI iteration also completed a timestep
   */
  virtual void _update_(bool timestepComplete, const Eigen::VectorXd& oldValues, const Eigen::VectorXd& res);


  static logging::Logger _log;

  // constant user-defined factors to scale the quasi-Newton system
  std::vector<double> _factors;

};

}}} // namespace precice, cplscheme, impl
