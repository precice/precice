#pragma once

#include <string>
#include "Action.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice::action {

class ScaleByAreaAction : public Action {
public:
  enum Scaling {
    /// Divides the data by the area of neighboring edges/triangles.
    SCALING_DIVIDE_BY_AREA,
    /// Multiplies the data by the area of neighboring edges/triangles.
    SCALING_MULTIPLY_BY_AREA
  };

  /**
   * @brief Constructor.
   *
   * @param[in] data Data that should be scaled.
   * @param[in] scalingType Type of scaling to be performed.
   */
  ScaleByAreaAction(
      Timing               timing,
      int                  targetDataID,
      const mesh::PtrMesh &mesh,
      Scaling              scaling);

  /**
   * @brief Scales data on mesh nodes according to selected scaling type.
   */
  void performAction() final;

private:
  logging::Logger _log{"action::ScaleByAreaAction"};

  mesh::PtrData _targetData;

  Scaling _scaling;
};

} // namespace precice::action
