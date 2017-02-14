#include "Geometry.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
#include "math/math.hpp"

namespace precice {
namespace geometry {

logging::Logger Geometry:: _log ( "precice::geometry::Geometry" );

Geometry:: Geometry
(
  const Eigen::VectorXd& offset )
:
  _offset (offset)
{}

void Geometry:: create
(
  mesh::Mesh& seed )
{
  TRACE(seed.getName());
  assertion ( seed.getDimensions() == _offset.size(), seed.getDimensions(),
               _offset.size() );
  specializedCreate ( seed );
  
  if ( not math::equals(getOffset(), Eigen::VectorXd::Zero(seed.getDimensions())) ) {
    Eigen::VectorXd temp( seed.getDimensions() );
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
