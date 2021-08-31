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
    double timeStepSize,
    double computedTimeWindowPart,
    double timeWindowSize)
{
  PRECICE_TRACE();
  const int       meshDimensions  = getMesh()->getDimensions();
  auto &          targetValues    = _targetData->values();
  const int       valueDimensions = _targetData->getDimensions();
  Eigen::VectorXd areas           = Eigen::VectorXd::Zero(getMesh()->vertices().size());
  PRECICE_ASSERT(targetValues.size() / valueDimensions == areas.size());

  if (meshDimensions == 2) {
    for (mesh::Edge &edge : getMesh()->edges()) {
      areas[edge.vertex(0).getID()] += edge.getEnclosingRadius();
      areas[edge.vertex(1).getID()] += edge.getEnclosingRadius();
    }
  } else {
    for (mesh::Triangle &face : getMesh()->triangles()) {
      areas[face.vertex(0).getID()] += face.getArea() / 3.0;
      areas[face.vertex(1).getID()] += face.getArea() / 3.0;
      areas[face.vertex(2).getID()] += face.getArea() / 3.0;
    }
  }
  if (_scaling == SCALING_DIVIDE_BY_AREA) {
    for (int i = 0; i < areas.size(); i++) {
      for (int dim = 0; dim < valueDimensions; dim++) {
        int valueIndex = i * valueDimensions + dim;
        targetValues[valueIndex] /= areas[i];
      }
    }
  } else if (_scaling == SCALING_MULTIPLY_BY_AREA) {
    for (int i = 0; i < areas.size(); i++) {
      for (int dim = 0; dim < valueDimensions; dim++) {
        int valueIndex = i * valueDimensions + dim;
        targetValues[valueIndex] *= areas[i];
      }
    }
  }
}

} // namespace action
} // namespace precice
