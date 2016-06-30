#ifndef PRECICE_ACTION_SCALEBYAREAACTION_HPP_
#define PRECICE_ACTION_SCALEBYAREAACTION_HPP_

#include "Action.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "tarch/logging/Log.h"

namespace precice {
  namespace mesh {
    class Edge;
    class Triangle;
  }
  namespace spacetree {
    class Spacetree;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace action {

class ScaleByAreaAction : public Action
{
public:

  enum Scaling {
    // @brief Divides the data by the area of neighboring edges/triangles.
    SCALING_DIVIDE_BY_AREA,
    // @brief Multiplies the data by the area of neighboring edges/triangles.
    SCALING_MULTIPLY_BY_AREA
  };

  /**
   * @brief Constructor.
   *
   * @param data [IN] Data that should be scaled.
   * @param scalingType [IN] Type of scaling to be performed.
   */
  ScaleByAreaAction (
    Timing               timing,
    int                  targetDataID,
    const mesh::PtrMesh& mesh,
    Scaling              scaling );

  /**
   * @brief Destructor.
   */
  virtual ~ScaleByAreaAction() {}

  /**
   * @brief Scales data on mesh nodes according to selected scaling type.
   *
   * At the moment, only a division of a property value by the associated area
   * of the neighboring edges (2D) is possible.
   */
  virtual void performAction (
    double time,
    double dt,
    double computedPartFullDt,
    double fullDt );

private:

  // @brief Logging device
  static tarch::logging::Log _log;

  mesh::PtrData _targetData;

  Scaling _scaling;
};

}} // namespace precice, action

#endif /* PRECICE_ACTION_SCALEBYAREAACTION_HPP_ */
