#pragma once

#include "action/Action.hpp"
#include "logging/Logger.hpp"

namespace precice {
  namespace tests {
    class BalanceVertexPositionActionTest;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace action {

/**
 * @brief iteratively balances the distance of vertices to each other.
 */
class BalanceVertexPositionAction : public Action
{
public:

  BalanceVertexPositionAction (
    Timing               timing,
    const mesh::PtrMesh& mesh,
    double               convergenceLimit,
    int                  maxIterations );

  virtual ~BalanceVertexPositionAction () {}

  virtual void performAction (
    double time,
    double dt,
    double computedPartFullDt,
    double fullDt );

private:

  static logging::Logger _log;

  double _eps;

  int _maxIterations;

  friend class tests::BalanceVertexPositionActionTest;
};

}} // namespace precice, action
