#pragma once

#include <string>
#include "action/Action.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace action {

/// Computes the curvature of a mesh geometry.
class ComputeCurvatureAction : public Action {
public:
  /// Constructor. Curvature values are stored in scalar data with given ID.
  ComputeCurvatureAction(
      Timing               timing,
      int                  dataID,
      const mesh::PtrMesh &mesh);

  /// Destructor, empty.
  virtual ~ComputeCurvatureAction() {}

  /// Computes the curvature of the mesh geometry.
  virtual void performAction(
      double time,
      double dt,
      double computedPartFullDt,
      double fullDt);

private:
  logging::Logger _log{"action::ComputeCurvatureAction"};

  mesh::PtrData _data;
};

} // namespace action
} // namespace precice
