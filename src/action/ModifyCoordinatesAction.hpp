#pragma once

#include "Action.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice
{
namespace action
{

/// Modifies a mesh's coordinates by using a coupling data set.
class ModifyCoordinatesAction : public Action
{
public:
  enum Mode {
    ADD_TO_COORDINATES_MODE,
    SUBTRACT_FROM_COORDINATES_MODE
  };

  ModifyCoordinatesAction(
      Timing               timing,
      int                  dataID,
      const mesh::PtrMesh &mesh,
      Mode                 mode);

  virtual ~ModifyCoordinatesAction(){};

  virtual void performAction(
      double time,
      double dt,
      double computedPartFullDt,
      double fullDt);

private:
  logging::Logger _log{"action::ModifyCoordinatesAction"};

  mesh::PtrData _data;

  Mode _mode;
};

} // namespace action
} // namespace precice
