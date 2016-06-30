#ifndef PRECICE_ACTION_BALANCEVERTEXPOSITIONACTION_HPP_
#define PRECICE_ACTION_BALANCEVERTEXPOSITIONACTION_HPP_

#include "action/Action.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"

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

  static tarch::logging::Log _log;

  double _eps;

  int _maxIterations;

  friend class tests::BalanceVertexPositionActionTest;
};

}} // namespace precice, action

#endif /* PRECICE_ACTION_BALANCEVERTEXPOSITIONACTION_HPP_ */
