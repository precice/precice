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
 * @brief Preconditioner that uses the values from the previous timestep to scale the quasi-Newton system.
 */
class ValuePreconditioner : public Preconditioner
{
public:

  ValuePreconditioner(std::vector<int> dimensions);
  /**
   * @brief Destructor, empty.
   */
  virtual ~ValuePreconditioner() {}

  /**
   * @brief Update the scaling after every FSI iteration.
   *
   * @param timestepComplete [IN] True if this FSI iteration also completed a timestep
   */
  virtual void update(bool timestepComplete, DataValues& oldValues, DataValues& res);

private:

  static tarch::logging::Log _log;

  bool _firstTimestep;
};

}}} // namespace precice, cplscheme, impl
