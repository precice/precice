#include "ScaleByDtAction.hpp"
#include <Eigen/Core>
#include <memory>
#include "action/Action.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace action {

ScaleByDtAction::ScaleByDtAction(
    Timing               timing,
    int                  sourceDataID,
    int                  targetDataID,
    const mesh::PtrMesh &mesh,
    Scaling              scaling)
    : Action(timing, mesh),
      _sourceData(mesh->data(sourceDataID)),
      _targetData(mesh->data(targetDataID)),
      _scaling(scaling)
{
  PRECICE_ASSERT(_sourceData->getDimensions() == _targetData->getDimensions(),
                 _sourceData->getDimensions(), _targetData->getDimensions());
}

void ScaleByDtAction::performAction(
    double time,
    double timeStepSize,
    double computedTimeWindowPart,
    double timeWindowSize)
{
  PRECICE_TRACE(timeStepSize, computedTimeWindowPart, timeWindowSize);
  auto &sourceValues = _sourceData->values();
  auto &targetValues = _targetData->values();
  PRECICE_ASSERT(sourceValues.size() == targetValues.size(),
                 sourceValues.size(), targetValues.size());
  if (_scaling == SCALING_BY_TIME_STEP_TO_TIME_WINDOW_RATIO) {
    double scaling = timeStepSize / timeWindowSize;
    PRECICE_DEBUG("Scale by ratio: time step size / time window size = {}", scaling);
    for (int i = 0; i < targetValues.size(); i++) {
      targetValues[i] = sourceValues[i] * scaling;
    }
  } else if (_scaling == SCALING_BY_TIME_WINDOW_SIZE) {
    PRECICE_DEBUG("Scale by dt {}", timeWindowSize);
    for (int i = 0; i < targetValues.size(); i++) {
      targetValues[i] = sourceValues[i] * timeWindowSize;
    }
  } else {
    PRECICE_ASSERT(_scaling == SCALING_BY_COMPUTED_TIME_WINDOW_PART_RATIO, _scaling);
    double scaling = computedTimeWindowPart / timeWindowSize;
    PRECICE_DEBUG("Scale by ratio: computed part of time window / time window size = {}", scaling);
    for (int i = 0; i < targetValues.size(); i++) {
      targetValues[i] = sourceValues[i] * scaling;
    }
  }
}

} // namespace action
} // namespace precice
