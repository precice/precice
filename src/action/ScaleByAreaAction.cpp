#include "ScaleByAreaAction.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace action {

logging::Logger ScaleByAreaAction:: _log ( "precice::action::ScaleByAreaAction" );

ScaleByAreaAction:: ScaleByAreaAction
(
  Timing               timing,
  int                  targetDataID,
  const mesh::PtrMesh& mesh,
  Scaling              scaling )
:
  Action(timing, mesh),
  _targetData(mesh->data (targetDataID)),
  _scaling(scaling)
{}

void ScaleByAreaAction:: performAction
(
  double time,
  double dt,
  double computedPartFullDt,
  double fullDt )
{
  preciceTrace("performAction()");
  preciceCheck(getMesh()->getDimensions() == 2, "performAction()",
                 "Not implemented for dimension != 2!");
  auto& targetValues = _targetData->values();
  utils::DynVector areas(getMesh()->vertices().size(), 0.0);
  for (mesh::Edge & edge : getMesh()->edges()){
    areas[edge.vertex(0).getID()] += edge.getEnclosingRadius();
    areas[edge.vertex(1).getID()] += edge.getEnclosingRadius();
  }
  int dimensions = _targetData->getDimensions();
  assertion(targetValues.size() / dimensions == areas.size());
  if (_scaling == SCALING_DIVIDE_BY_AREA){
    for (int i=0; i < areas.size(); i++){
      for (int dim=0; dim < dimensions; dim++){
        int valueIndex = i*dimensions + dim;
        targetValues[valueIndex] /= areas[i];
      }
    }
  }
  else if (_scaling == SCALING_MULTIPLY_BY_AREA){
    for (int i=0; i < areas.size(); i++){
      for (int dim=0; dim < dimensions; dim++){
        int valueIndex = i*dimensions + dim;
        targetValues[valueIndex] *= areas[i];
      }
    }
  }
}

}} // namspace precice, action

