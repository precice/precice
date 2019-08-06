#include "ModifyCoordinatesAction.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"

namespace precice
{
namespace action
{

ModifyCoordinatesAction::ModifyCoordinatesAction(
    Timing               timing,
    int                  dataID,
    const mesh::PtrMesh &mesh,
    Mode                 mode)
    : Action(timing, mesh),
      _data(mesh->data(dataID)),
      _mode(mode)
{
}

void ModifyCoordinatesAction::performAction(
    double time,
    double dt,
    double computedPartFullDt,
    double fullDt)
{
  P_TRACE();
  auto &          values = _data->values();
  int             dim    = getMesh()->getDimensions();
  Eigen::VectorXd data(dim);
  if (_mode == ADD_TO_COORDINATES_MODE) {
    P_DEBUG("Adding data to coordinates");
    for (mesh::Vertex &vertex : getMesh()->vertices()) {
      for (int i = 0; i < dim; i++) {
        data[i] = values[vertex.getID() * dim + i];
      }
      data += vertex.getCoords();
      vertex.setCoords(data);
    }
  } else if (_mode == SUBTRACT_FROM_COORDINATES_MODE) {
    P_DEBUG("Subtracting data from coordinates");
    for (mesh::Vertex &vertex : getMesh()->vertices()) {
      data = vertex.getCoords();
      for (int i = 0; i < dim; i++) {
        data[i] -= values[vertex.getID() * dim + i];
      }
      vertex.setCoords(data);
    }
  } else {
    P_ERROR("Unknown mode type!");
  }
  getMesh()->computeState();
  getMesh()->meshChanged(*getMesh());
}

} // namespace action
} // namespace precice
