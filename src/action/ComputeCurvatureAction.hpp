// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_ACTION_COMPUTECURVATUREACTION_HPP_
#define PRECICE_ACTION_COMPUTECURVATUREACTION_HPP_

#include "action/Action.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"

namespace precice {
namespace action {

/**
 * @brief Computes the curvature of a mesh geometry.
 */
class ComputeCurvatureAction : public Action
{
public:

  /**
   * @brief Constructor. Curvature values are stored in scalar data with given ID.
   */
  ComputeCurvatureAction (
    Timing               timing,
    int                  dataID,
    const mesh::PtrMesh& mesh );

  /**
   * @brief Destructor, empty.
   */
  virtual ~ComputeCurvatureAction() {}

  /**
   * @brief Computes the curvature of the mesh geometry.
   */
  virtual void performAction (
    double time,
    double dt,
    double computedPartFullDt,
    double fullDt );

private:

  // @brief Logging device
  static tarch::logging::Log _log;

  mesh::PtrData _data;
};

}} // namespace precice, action

#endif /* PRECICE_ACTION_COMPUTECURVATUREACTION_HPP_ */
