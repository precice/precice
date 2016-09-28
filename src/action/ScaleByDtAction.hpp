#ifndef PRECICE_ACTION_SCALEBYDTACTION_HPP_
#define PRECICE_ACTION_SCALEBYDTACTION_HPP_

#include "Action.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "logging/Logger.hpp"

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

class ScaleByDtAction : public Action
{
public:

  enum Scaling {
    // @brief Scales data by ratio of last computed timestep to full timestep length.
    SCALING_BY_COMPUTED_DT_RATIO,
    // @brief Scales data by last computed timestep
    SCALING_BY_DT,
    // @brief Scales data by ratio of computed part of full timestep.
    SCALING_BY_COMPUTED_DT_PART_RATIO
  };

  /**
   * @brief Constructor.
   *
   * @param data [IN] Data that should be scaled.
   * @param scalingType [IN] Type of scaling to be performed.
   */
  ScaleByDtAction (
    Timing               timing,
    int                  sourceDataID,
    int                  targetDataID,
    const mesh::PtrMesh& mesh,
    Scaling              scaling );

  /**
   * @brief Destructor.
   */
  virtual ~ScaleByDtAction() {}

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
  static logging::Logger _log;

  mesh::PtrData _sourceData;

  mesh::PtrData _targetData;

  Scaling _scaling;
};

}} // namespace precice, action

#endif /* PRECICE_ACTION_SCALEBYDTACTION_HPP_ */
