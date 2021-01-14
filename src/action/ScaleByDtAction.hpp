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
    /// Scales data by ratio of last computed time step to time window size.
    SCALING_BY_TIME_STEP_TO_TIME_WINDOW_RATIO,
    /// Scales data by time window size.
    SCALING_BY_TIME_WINDOW_SIZE,
    /// Scales data by ratio of computed part of time window to length of time window.
    SCALING_BY_COMPUTED_TIME_WINDOW_PART_RATIO
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
   */
  virtual void performAction(
      double time,
      double timeStepSize,
      double computedTimeWindowPart,
      double timeWindowSize);

private:
  logging::Logger _log{"action::ScaleByDtAction"};

  mesh::PtrData _sourceData;

  mesh::PtrData _targetData;

  Scaling _scaling;
};

} // namespace action
} // namespace precice
