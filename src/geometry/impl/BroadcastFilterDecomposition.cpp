#include "BroadcastFilterDecomposition.hpp"
#include "com/Communication.hpp"
#include "com/CommunicateMesh.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"
#include "utils/EventTimings.hpp"

using precice::utils::Event;

namespace precice {
namespace geometry {
namespace impl {

tarch::logging::Log BroadcastFilterDecomposition:: _log ( "precice::geometry::BroadcastFilterDecomposition" );

BroadcastFilterDecomposition:: BroadcastFilterDecomposition
(
  int    dimensions)
  :
  Decomposition ( dimensions )
{}


void BroadcastFilterDecomposition:: decompose(
  mesh::Mesh& seed)
{
  preciceTrace1 ( "decompose()", utils::MasterSlave::_rank );
  using tarch::la::raw;

  std::map<int,std::vector<int> > boundingVertexDistribution;
  std::vector<int> filteredVertexPositions;

  broadcast(seed);
  filter(seed, filteredVertexPositions);
  feedback(seed, filteredVertexPositions);
}

void BroadcastFilterDecomposition:: broadcast(
  mesh::Mesh& seed)
{
  preciceTrace1 ( "broadcast()", utils::MasterSlave::_rank );
  preciceInfo("broadcast()", "Broadcast mesh " << seed.getName() );
  Event e("broadcast mesh");


  if (utils::MasterSlave::_slaveMode) {
    com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh (seed, 0);
  }
  else{ // Master
    assertion(utils::MasterSlave::_rank==0);
    assertion(utils::MasterSlave::_size>1);
    seed.setGlobalNumberOfVertices(seed.vertices().size());

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
      com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh ( seed, rankSlave );
    }
  }
}

void BroadcastFilterDecomposition:: filter(
  mesh::Mesh& seed,
  std::vector<int>& filteredVertexPositions)
{
  preciceTrace1 ( "filter()", utils::MasterSlave::_rank );
  preciceInfo("filter()", "Filter mesh " << seed.getName() );
  Event e("filter mesh");

  seed.computeState();
  computeBoundingMappings();
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());
  filteredVertexPositions = filterMesh(seed, filteredMesh);
  seed.clear();
  seed.addMesh(filteredMesh);
  clearBoundingMappings();
}

void BroadcastFilterDecomposition:: feedback(
  mesh::Mesh& seed,
  std::vector<int>& filteredVertexPositions)
{
  preciceTrace1 ( "feedback()", utils::MasterSlave::_rank );
  preciceInfo("feedback()", "Feedback mesh " << seed.getName() );
  Event e("feedback mesh");

  int numberOfVertices = filteredVertexPositions.size();

  if (utils::MasterSlave::_slaveMode) {
    utils::MasterSlave::_communication->send(numberOfVertices,0);
    if (numberOfVertices!=0) {
      utils::MasterSlave::_communication->send(filteredVertexPositions.data(),numberOfVertices,0);
    }
  }
  else { // Master
    assertion(utils::MasterSlave::_rank==0);
    assertion(utils::MasterSlave::_size>1);

    seed.getVertexDistribution()[0] = filteredVertexPositions;

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      int numberOfVertices = -1;
      utils::MasterSlave::_communication->receive(numberOfVertices,rankSlave);
      std::vector<int> slaveVertexIDs(numberOfVertices,-1);
      if (numberOfVertices!=0) {
        utils::MasterSlave::_communication->receive(slaveVertexIDs.data(),numberOfVertices,rankSlave);
      }

      seed.getVertexDistribution()[rankSlave] = slaveVertexIDs;
    }
  }

}

bool BroadcastFilterDecomposition:: doesVertexContribute(
  const mesh::Vertex& vertex)
{
  //works as easy as this since only read-consistent and write-conservative are allowed
  assertion(_boundingFromMapping.use_count()>0 || _boundingToMapping.use_count()>0);
  bool exit = false;
  if (_boundingFromMapping.use_count() > 0) {
    exit = exit || _boundingFromMapping->doesVertexContribute(vertex.getID());
  }
  if (_boundingToMapping.use_count() > 0) {
    exit = exit || _boundingToMapping->doesVertexContribute(vertex.getID());
  }
  return exit;
}

}}} // namespace precice, geometry, impl
