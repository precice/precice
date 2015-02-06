// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SolverGeometry.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace geometry {

tarch::logging::Log SolverGeometry:: _log ( "precice::geometry::SolverGeometry" );

SolverGeometry:: SolverGeometry
(
  const utils::DynVector&  offset)
:
  Geometry ( offset )
{
  preciceTrace ( "SolverGeometry()" );
}


void SolverGeometry:: specializedCreate
(
  mesh::Mesh& seed )
{
  preciceTrace1 ( "specializedCreate()", seed.getName() );
  //nothing to do here
}

}} // namespace precice, geometry
