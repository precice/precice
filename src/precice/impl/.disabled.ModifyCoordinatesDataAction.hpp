// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_MODIFYCOORDINATESDATAACTION_HPP_
#define PRECICE_IMPL_MODIFYCOORDINATESDATAACTION_HPP_

#include "AbstractDataAction.hpp"
#include "MeshContext.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"

namespace precice {
namespace impl {

/**
 * @brief Modifies a mesh's coordinates by using a coupling data set.
 */
class ModifyCoordinatesDataAction : public AbstractDataAction
{
public:

  enum ModeConstants {
    ADD_TO_COORDINATES_MODE,
    SUBTRACT_FROM_COORDINATES_MODE
  };

  ModifyCoordinatesDataAction (
    TimingConstants timing,
    int             dataID,
    MeshContext &   meshContext,
    ModeConstants   mode );

  virtual ~ModifyCoordinatesDataAction() {};

  virtual void performAction (
    double dt,
    double computedPartFullDt,
    double fullDt );

private:

  static tarch::logging::Log _log;

  mesh::PtrData _data;

  MeshContext & _meshContext;

  ModeConstants _mode;
};


}} // namespace precice, impl

#endif /* PRECICE_IMPL_MODIFYCOORDINATESDATAACTION_HPP_ */
