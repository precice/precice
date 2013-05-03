// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_COMPUTECURVATUREDATAACTION_HPP_
#define PRECICE_IMPL_COMPUTECURVATUREDATAACTION_HPP_

#include "precice/impl/AbstractDataAction.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace impl {

/**
 * @brief Computes the curvature of a mesh geometry.
 */
class ComputeCurvatureDataAction : public AbstractDataAction
{
public:

  /**
   * @brief Constructor. Curvature values are stored in scalar data with given ID.
   */
  ComputeCurvatureDataAction (
    TimingConstants       timing,
    int                   dataID,
    const mesh::PtrMesh & mesh );

  /**
   * @brief Destructor, empty.
   */
  virtual ~ComputeCurvatureDataAction () {}

  /**
   * @brief Computes the curvature of the mesh geometry.
   */
  virtual void performAction (
    double dt,
    double computedPartFullDt,
    double fullDt );

private:

  // @brief Logging device
  static tarch::logging::Log _log;

  mesh::PtrData _data;

  mesh::PtrMesh _mesh;
};

}} // namespace precice, impl

#endif /* PRECICE_IMPL_COMPUTECURVATUREDATAACTION_HPP_ */
