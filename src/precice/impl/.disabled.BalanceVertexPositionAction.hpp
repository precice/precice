// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_BALANCEVERTEXPOSITIONACTION_HPP_
#define PRECICE_IMPL_BALANCEVERTEXPOSITIONACTION_HPP_

#include "precice/impl/AbstractDataAction.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"

namespace precice {
  namespace tests {
    class BalanceVertexPositionActionTest;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace impl {

/**
 * @brief iteratively balances the distance of vertices to each other.
 */
class BalanceVertexPositionAction : public AbstractDataAction
{
public:

  BalanceVertexPositionAction (
    TimingConstants       timing,
    const mesh::PtrMesh & mesh,
    double                convergenceLimit,
    int                   maxIterations );

  virtual ~BalanceVertexPositionAction () {}

  virtual void performAction (
    double dt,
    double computedPartFullDt,
    double fullDt );

private:

  static tarch::logging::Log _log;

  mesh::PtrMesh _mesh;

  double _eps;

  int _maxIterations;

  friend class tests::BalanceVertexPositionActionTest;
};

}} // namespace precice, impl

#endif /* PRECICE_IMPL_BALANCEVERTEXPOSITIONACTION_HPP_ */
