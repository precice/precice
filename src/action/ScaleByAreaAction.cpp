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

namespace precice::action {

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

void ScaleByAreaAction::performAction()
{
  PRECICE_TRACE();
  const int meshDimensions = getMesh()->getDimensions();

  for (auto &targetStample : _targetData->stamples()) {

    auto &targetValues              = _targetData->values();
    targetValues                    = targetStample.sample.values;
    const int       valueDimensions = _targetData->getDimensions();
    Eigen::VectorXd areas           = Eigen::VectorXd::Zero(getMesh()->vertices().size());
    PRECICE_ASSERT(targetValues.size() / valueDimensions == areas.size());

    if (meshDimensions == 2) {
      PRECICE_CHECK(getMesh()->edges().size() != 0,
                    "The multiply/divide-by-area actions require meshes with connectivity information. In 2D, please ensure that the mesh {} contains edges.", getMesh()->getName());
      for (mesh::Edge &edge : getMesh()->edges()) {
        areas[edge.vertex(0).getID()] += edge.getEnclosingRadius();
        areas[edge.vertex(1).getID()] += edge.getEnclosingRadius();
      }
    } else {
      PRECICE_CHECK(getMesh()->triangles().size() != 0,
                    "The multiply/divide-by-area actions require meshes with connectivity information. In 3D, please ensure that the mesh {} contains triangles.", getMesh()->getName());
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
    _targetData->setSampleAtTime(targetStample.timestamp, _targetData->sample());
  }
}

} // namespace precice::action
