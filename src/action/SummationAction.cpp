#include "SummationAction.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "action/Action.hpp"
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"

namespace precice::action {

SummationAction::SummationAction(
    Timing                  timing,
    const std::vector<int> &sourceDataIDs,
    int                     targetDataID,
    const mesh::PtrMesh &   mesh)
    : Action(timing, mesh, mapping::Mapping::MeshRequirement::VERTEX), _targetData(mesh->data(targetDataID))
{

  for (int sourceID : sourceDataIDs) {
    _sourceDataVector.push_back(mesh->data(sourceID));
  }

  for (const auto &source : _sourceDataVector) {
    PRECICE_CHECK(source->getDimensions() == _targetData->getDimensions(), "Source and target data dimensions (scalar or vector) of summation action need to be identical.");
  }
}

void SummationAction::performAction()
{
  PRECICE_TRACE();

  const auto &referenceData = _sourceDataVector.front(); // serves as reference
  const int   nStamples     = referenceData->stamples().size();
  for (int stampleId = 0; stampleId < nStamples; stampleId++) { // simultaneously loop over stamples in all sourceData of _sourceDataVector
    auto &targetValues = _targetData->values();
    targetValues.setZero();
    const double currentTimestamp = referenceData->stamples()[stampleId].timestamp;
    for (const auto &sourceData : _sourceDataVector) {
      auto sourceStample = sourceData->stamples()[stampleId];
      PRECICE_CHECK(math::equals(sourceStample.timestamp, currentTimestamp), "Trying to perform summation action on samples with different timestamps: expected timestamp {}, but got source data with timestamp {}.  Time meshes of all source data must agree. Actions do not fully support subcycling yet.", currentTimestamp, sourceStample.timestamp);
      auto sourceDataValues = sourceStample.sample.values;
      targetValues += sourceDataValues;
    }
    _targetData->setSampleAtTime(currentTimestamp, _targetData->sample());
  }
}

} // namespace precice::action
