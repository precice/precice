#include "ScaleByAreaAction.hpp"
#include <Eigen/Core>
#include <memory>
#include "action/Action.hpp"
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace action {

ScaleByAreaAction::ScaleByAreaAction(
    Timing               timing,
    int                  targetDataID,
    const mesh::PtrMesh &mesh,
    Scaling              scaling)
    : Action(timing, mesh, mapping::Mapping::MeshRequirement::FULL),
      _targetData(mesh->data(targetDataID)),
      _scaling(scaling)
{
}

void ScaleByAreaAction::performAction(
    double time,
    double dt,
    double computedPartFullDt,
    double fullDt)
{
  PRECICE_TRACE();
  PRECICE_ASSERT(getMesh()->getDimensions() == 2, "The ScaleByAreaAction requires a mesh of dimensionality 2.");
  auto &          targetValues = _targetData->values();
  Eigen::VectorXd areas        = Eigen::VectorXd::Zero(getMesh()->vertices().size());
  for (mesh::Edge &edge : getMesh()->edges()) {
    areas[edge.vertex(0).getID()] += edge.getEnclosingRadius();
    areas[edge.vertex(1).getID()] += edge.getEnclosingRadius();
  }
  int dimensions = _targetData->getDimensions();
  PRECICE_ASSERT(targetValues.size() / dimensions == areas.size());
  if (_scaling == SCALING_DIVIDE_BY_AREA) {
    for (int i = 0; i < areas.size(); i++) {
      for (int dim = 0; dim < dimensions; dim++) {
        int valueIndex = i * dimensions + dim;
        targetValues[valueIndex] /= areas[i];
      }
    }
  } else if (_scaling == SCALING_MULTIPLY_BY_AREA) {
    for (int i = 0; i < areas.size(); i++) {
      for (int dim = 0; dim < dimensions; dim++) {
        int valueIndex = i * dimensions + dim;
        targetValues[valueIndex] *= areas[i];
      }
    }
  }
}

} // namespace action
} // namespace precice
