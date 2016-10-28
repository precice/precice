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

logging::Logger BroadcastFilterDecomposition:: _log ( "precice::geometry::BroadcastFilterDecomposition" );

BroadcastFilterDecomposition:: BroadcastFilterDecomposition
(
  int    dimensions, double safetyFactor )
  :
  Decomposition ( dimensions, safetyFactor  )
{}


void BroadcastFilterDecomposition:: decompose(
  mesh::Mesh& seed)
{
  TRACE(utils::MasterSlave::_rank );

  std::map<int,std::vector<int> > boundingVertexDistribution;
  std::vector<int> filteredVertexPositions;

  broadcast(seed);
  filter(seed, filteredVertexPositions);
  feedback(seed, filteredVertexPositions);
}

void BroadcastFilterDecomposition:: broadcast(
  mesh::Mesh& seed)
{
  preciceTrace ( "broadcast()", utils::MasterSlave::_rank );
  preciceInfo("broadcast()", "Broadcast mesh " << seed.getName() );
  Event e("broadcast mesh");


  if (utils::MasterSlave::_slaveMode) {
    com::CommunicateMesh(utils::MasterSlave::_communication).broadcastReceiveMesh (seed);
  }
  else{ // Master
    assertion(utils::MasterSlave::_rank==0);
    assertion(utils::MasterSlave::_size>1);
    com::CommunicateMesh(utils::MasterSlave::_communication).broadcastSendMesh ( seed);
  }
}

void BroadcastFilterDecomposition:: filter(
  mesh::Mesh& seed,
  std::vector<int>& filteredVertexPositions)
{
  preciceTrace ( "filter()", utils::MasterSlave::_rank );
  preciceInfo("filter()", "Filter mesh " << seed.getName() );
  Event e("filter mesh");

  // first, bounding box filter
  DEBUG("First Filter BB, #vertices " << seed.vertices().size());
  assertion(not _filterByMapping);
  _bb = mesh::Mesh::BoundingBox (_dimensions,
                   std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));
  mergeBoundingBoxes(_bb);
  for (int d=0; d<_dimensions; d++) {
    if (_bb[d].second > _bb[d].first && _safetyGap < _bb[d].second - _bb[d].first)
      _safetyGap = _bb[d].second - _bb[d].first;
  }
  assertion(_safetyFactor>=0.0);
  _safetyGap *= _safetyFactor;
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());
  std::vector<int> tmpVertexPostitions = filterMesh(seed, filteredMesh);

  // second, mapping filter
  DEBUG("Second Filter Mapping, #vertices " << filteredMesh.vertices().size());
  _filterByMapping = true;
  seed.clear();
  seed.addMesh(filteredMesh);
  seed.computeState();
  computeBoundingMappings();
  filteredMesh.clear();
  filteredVertexPositions = filterMesh(seed, filteredMesh);
  seed.clear();
  seed.addMesh(filteredMesh);
  clearBoundingMappings();

  //merge the 2 filters
  DEBUG("Merge Filters, #vertices " << filteredMesh.vertices().size());
  for(size_t i=0;i<filteredVertexPositions.size();i++){
    filteredVertexPositions[i] = tmpVertexPostitions[filteredVertexPositions[i]];
  }

}

void BroadcastFilterDecomposition:: feedback(
  mesh::Mesh& seed,
  std::vector<int>& filteredVertexPositions)
{
  preciceTrace ( "feedback()", utils::MasterSlave::_rank );
  preciceInfo("feedback()", "Feedback mesh " << seed.getName() );

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


}}} // namespace precice, geometry, impl
