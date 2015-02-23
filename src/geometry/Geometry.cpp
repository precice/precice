// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Geometry.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/Globals.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace geometry {

tarch::logging::Log Geometry:: _log ( "precice::geometry::Geometry" );

Geometry:: Geometry
(
  const utils::DynVector& offset )
:
  _offset (offset)
{}

void Geometry:: create
(
  mesh::Mesh& seed )
{
  preciceTrace1 ( "create()", seed.getName() );
  assertion2 ( seed.getDimensions() == _offset.size(), seed.getDimensions(),
               _offset.size() );
  specializedCreate ( seed );
  utils::DynVector zero ( seed.getDimensions(), 0.0 );
  if ( not tarch::la::equals(getOffset(), zero) ) {
    utils::DynVector temp ( seed.getDimensions() );
    for (mesh::Vertex& vertex : seed.vertices()) {
      temp = _offset;
      temp += vertex.getCoords();
      vertex.setCoords ( temp );
    }
  }
  seed.computeState();
  allocateDataValues(seed);
}

void Geometry:: allocateDataValues ( mesh::Mesh & mesh )
{
  mesh.allocateDataValues ();
}

}} // namespace precice, geometry
