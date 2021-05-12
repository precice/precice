#include "SummationAction.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "action/Action.hpp"
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace action {

SummationAction::SummationAction(
    Timing               timing,
    std::vector<int>     sourceDataIDs,
    int                  targetDataID,
    const mesh::PtrMesh &mesh)
    : Action(timing, mesh, mapping::Mapping::MeshRequirement::VERTEX), _targetData(mesh->data(targetDataID))
{

  for (int sourceID : sourceDataIDs) {
    _sourceDataVector.push_back(mesh->data(sourceID));
  }

  for (const auto &source : _sourceDataVector) {
    PRECICE_CHECK(source->getDimensions() == _targetData->getDimensions(), "Source and target data dimensions (scalar or vector) of summation action need to be identical.");
  }
}

void SummationAction::performAction(
    double time,
    double timeStepSize,
    double computedTimeWindowPart,
    double timeWindowSize)
{
  PRECICE_TRACE();

  auto &targetValues = _targetData->values();
  targetValues.setZero();

  for (const auto &sourceData : _sourceDataVector) {
    targetValues += sourceData->values();
  }
}

} // namespace action
} // namespace precice