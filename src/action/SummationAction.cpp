#include "SummationAction.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <memory>
#include "action/Action.hpp"
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "utils/ranges.hpp"

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

  // check if all source stamples have identical times
  const auto &referenceData = _sourceDataVector.front(); // serves as reference
  auto        times         = referenceData->timeStepsStorage().getTimes();
  for (const auto &sourceData : utils::ranges::tail(_sourceDataVector)) {
    auto sourceTimes = sourceData->timeStepsStorage().getTimes();
    PRECICE_CHECK(times == sourceTimes,
                  "Trying to perform summation action on data with different timestamps data \"{}\" has {}, while data \"{}\" has {}. Time meshes of all source data must agree. Actions do not fully support subcycling yet.",
                  referenceData->getName(), times, sourceData->getName(), sourceTimes);
  }

  // Sum up the samples, note that _targetData does not contain stamples
  for (int stampleId = 0, nStamples = referenceData->stamples().size(); stampleId < nStamples; stampleId++) { // simultaneously loop over stamples in all sourceData of _sourceDataVector
    time::Stample sum = _sourceDataVector.front()->stamples()[stampleId];
    for (const auto &sourceData : utils::ranges::tail(_sourceDataVector)) {
      sum.sample.values += sourceData->stamples()[stampleId].sample.values;
    }
    _targetData->setSampleAtTime(sum.timestamp, std::move(sum.sample));
  }
}

} // namespace precice::action
