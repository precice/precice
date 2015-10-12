// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_GEOMETRY_SOLVERGEOMETRY_HPP_
#define PRECICE_GEOMETRY_SOLVERGEOMETRY_HPP_

#include "Geometry.hpp"
#include "com/Communication.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"
#include <string>
#include <map>

namespace precice {
namespace geometry {

/**
 * @brief Geometry for those meshes that are not communicated (meshes with "provide"
 * but no counterpart "from"). This class is a simple wrapper without any functionality.
 */
class SolverGeometry : public Geometry
{
public:

  SolverGeometry (
    const utils::DynVector& offset);

  virtual ~SolverGeometry() {}

protected:

  /**
   * @brief Is called from Geometry. Nothing to be done here.
   */
  virtual void specializedCreate ( mesh::Mesh& seed );

private:

  // @brief Logging device.
  static tarch::logging::Log _log;
};

}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_SOLVERGEOMETRY_HPP_ */
