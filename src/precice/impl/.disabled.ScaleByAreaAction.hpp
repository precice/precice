// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_SCALEBYAREAACTION_HPP_
#define PRECICE_IMPL_SCALEBYAREAACTION_HPP_

#include "AbstractDataAction.hpp"
#include "MeshContext.hpp"
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
namespace impl {

class ScaleByAreaAction : public AbstractDataAction
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
    TimingConstants       timing,
    int                   targetDataID,
    const mesh::PtrMesh & mesh,
    Scaling               scaling );

  /**
   * @brief Destructor.
   */
  virtual ~ScaleByAreaAction () {}

  /**
   * @brief Scales data on mesh nodes according to selected scaling type.
   *
   * At the moment, only a division of a property value by the associated area
   * of the neighboring edges (2D) is possible.
   */
  virtual void performAction (
    double dt,
    double computedPartFullDt,
    double fullDt );

private:

  // @brief Logging device
  static tarch::logging::Log _log;

  mesh::PtrData _targetData;

  mesh::PtrMesh _mesh;

  Scaling _scaling;
};

}} // namespace precice, impl

#endif /* PRECICE_IMPL_SCALEBYAREAACTION_HPP_ */
