#include "SummationAction.hpp"

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
    PRECICE_CHECK(source->getDimensions() == _targetData->getDimensions(), "Source and target data dimensions should be same for summation action.");
  }
}

void SummationAction::performAction(
    double time,
    double dt,
    double computedPartFullDt,
    double fullDt)
{
  PRECICE_TRACE();
  auto &targetValues = _targetData->values();
  auto  sum          = _sourceDataVector.at(0)->values();

  for (int i = 1; i < _sourceDataVector.size(); ++i) {
    sum += _sourceDataVector.at(i)->values();
  }
  targetValues = sum;
}

} // namespace action
} // namespace precice