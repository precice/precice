#pragma once

#include <string>
#include "Action.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace action {

class ScaleByDtAction : public Action {
public:
  enum Scaling {
    /// Scales data by ratio of last computed timestep to full timestep length.
    SCALING_BY_COMPUTED_DT_RATIO,
    /// Scales data by last computed timestep
    SCALING_BY_DT,
    /// Scales data by ratio of computed part of full timestep.
    SCALING_BY_COMPUTED_DT_PART_RATIO
  };

  /**
   * @brief Constructor.
   *
   * @param[in] data Data that should be scaled.
   * @param[in] scalingType Type of scaling to be performed.
   */
  ScaleByDtAction(
      Timing               timing,
      int                  sourceDataID,
      int                  targetDataID,
      const mesh::PtrMesh &mesh,
      Scaling              scaling);

  virtual ~ScaleByDtAction() {}

  /**
   * @brief Scales data on mesh nodes according to selected scaling type.
   *
   * At the moment, only a division of a property value by the associated area
   * of the neighboring edges (2D) is possible.
   */
  virtual void performAction(
      double time,
      double dt,
      double computedPartFullDt,
      double fullDt);

private:
  logging::Logger _log{"action::ScaleByDtAction"};

  mesh::PtrData _sourceData;

  mesh::PtrData _targetData;

  Scaling _scaling;
};

} // namespace action
} // namespace precice
