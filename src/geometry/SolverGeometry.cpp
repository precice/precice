#include "SolverGeometry.hpp"
#include "com/Communication.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace geometry {

logging::Logger SolverGeometry:: _log ( "precice::geometry::SolverGeometry" );

SolverGeometry:: SolverGeometry
(
  const Eigen::VectorXd& offset)
:
  Geometry ( offset )
{
  TRACE();
}


void SolverGeometry:: specializedCreate
(
  mesh::Mesh& seed )
{
  TRACE(seed.getName());

  //generate vertexDistribution also for non-communicated geometries as this information
  //is needed to assign global indices resp. vertexOffsets

  int numberOfVertices = seed.vertices().size();

  if (utils::MasterSlave::_slaveMode) {
    utils::MasterSlave::_communication->send(numberOfVertices,0);
  }
  else if(utils::MasterSlave::_masterMode) {
    assertion(utils::MasterSlave::_size>1);
    int vertexCounter = 0;

    //add master vertices
    for(int i=0; i<numberOfVertices; i++){
      seed.getVertexDistribution()[0].push_back(vertexCounter);
      vertexCounter++;
    }

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      utils::MasterSlave::_communication->receive(numberOfVertices,rankSlave);
      for(int i=0; i<numberOfVertices; i++){
        seed.getVertexDistribution()[rankSlave].push_back(vertexCounter);
        vertexCounter++;
      }
    }
    seed.setGlobalNumberOfVertices(vertexCounter);
  }
  else{ //coupling mode
    seed.setGlobalNumberOfVertices(seed.vertices().size());
  }
}

}} // namespace precice, geometry
