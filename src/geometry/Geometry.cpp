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
  assertion ( seed.getDimensions() == _offset.size(), seed.getDimensions(),
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
