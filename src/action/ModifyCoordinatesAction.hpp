// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_ACTION_MODIFYCOORDINATESACTION_HPP_
#define PRECICE_ACTION_MODIFYCOORDINATESACTION_HPP_

#include "Action.hpp"
#include "mesh/SharedPointer.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace action {

/**
 * @brief Modifies a mesh's coordinates by using a coupling data set.
 */
class ModifyCoordinatesAction : public Action
{
public:

  enum Mode {
    ADD_TO_COORDINATES_MODE,
    SUBTRACT_FROM_COORDINATES_MODE
  };

  ModifyCoordinatesAction (
    Timing               timing,
    int                  dataID,
    const mesh::PtrMesh& mesh,
    Mode                 mode );

  virtual ~ModifyCoordinatesAction() {};

  virtual void performAction (
    double time,
    double dt,
    double computedPartFullDt,
    double fullDt );

private:

  static logging::Logger _log;

  mesh::PtrData _data;

  Mode _mode;
};


}} // namespace precice, action

#endif /* PRECICE_ACTION_MODIFYCOORDINATESACTION_HPP_ */
