#include "PreFilterPostFilterDecomposition.hpp"
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

tarch::logging::Log PreFilterPostFilterDecomposition:: _log ( "precice::geometry::PreFilterPostFilterDecomposition" );

PreFilterPostFilterDecomposition:: PreFilterPostFilterDecomposition
(
  int    dimensions,
  double safetyFactor)
  :
  Decomposition ( dimensions, safetyFactor )
{}


void PreFilterPostFilterDecomposition:: decompose(
  mesh::Mesh& seed)
{
  preciceTrace1 ( "decompose()", utils::MasterSlave::_rank );
  using tarch::la::raw;

  std::map<int,std::vector<int> > boundingVertexDistribution;
  std::vector<int> filteredVertexPositions;

  preFilter(seed, boundingVertexDistribution);
  _filterByMapping = true;
  postFilter(seed, filteredVertexPositions);
  feedback(seed, boundingVertexDistribution, filteredVertexPositions);
}

void PreFilterPostFilterDecomposition:: preFilter(
  mesh::Mesh& seed,
  std::map<int,std::vector<int> >& boundingVertexDistribution)
{
  preciceTrace1 ( "preFilter()", utils::MasterSlave::_rank );
  preciceInfo("preFilter()", "Pre-filter mesh " << seed.getName() );
  Event e("pre-filter mesh");

  assertion(not _filterByMapping);

  if (utils::MasterSlave::_slaveMode) {
    mesh::Mesh::BoundingBox bb = mesh::Mesh::BoundingBox (_dimensions,
                             std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));
    mergeBoundingBoxes(bb);
    com::CommunicateMesh(utils::MasterSlave::_communication).sendBoundingBox (bb, 0);
    com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh (seed, 0);
  }
  else{ // Master
    assertion(utils::MasterSlave::_rank==0);
    assertion(utils::MasterSlave::_size>1);

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
      _bb = mesh::Mesh::BoundingBox (_dimensions, std::make_pair(0.0,0.0));
      com::CommunicateMesh(utils::MasterSlave::_communication).receiveBoundingBox ( _bb, rankSlave);

      for (int d=0; d<_dimensions; d++) {
        if (_bb[d].second > _bb[d].first && _safetyGap < _bb[d].second - _bb[d].first)
          _safetyGap = _bb[d].second - _bb[d].first;
      }
      assertion(_safetyFactor>=0.0);
      _safetyGap *= _safetyFactor;

      preciceDebug("From slave " << rankSlave << ", bounding mesh: " << _bb[0].first
                   << ", " << _bb[0].second << " and " << _bb[1].first << ", " << _bb[1].second);
      mesh::Mesh slaveMesh("SlaveMesh", _dimensions, seed.isFlipNormals());
      boundingVertexDistribution[rankSlave] = filterMesh(seed, slaveMesh);
      com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh ( slaveMesh, rankSlave );
    }

    // Now also filter the remaining master mesh
    _bb = mesh::Mesh::BoundingBox (_dimensions,
               std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));
    mergeBoundingBoxes(_bb);
    for (int d=0; d < _dimensions; d++) {
      if (_bb[d].second > _bb[d].first && _safetyGap < _bb[d].second - _bb[d].first)
        _safetyGap = _bb[d].second - _bb[d].first;
    }
    _safetyGap *= _safetyFactor;
    mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());
    boundingVertexDistribution[0] = filterMesh(seed, filteredMesh);
    seed.clear(); //clear global mesh on master
    seed.addMesh(filteredMesh);
    preciceDebug("Master mesh after filtering, #vertices " << seed.vertices().size());
  }
}

void PreFilterPostFilterDecomposition:: postFilter(
  mesh::Mesh& seed,
  std::vector<int>& filteredVertexPositions)
{
  preciceTrace1 ( "postFilter()", utils::MasterSlave::_rank );
  preciceInfo("postFilter()", "Post-filter mesh " << seed.getName() );
  Event e("post-filter mesh");

  seed.computeState();
  computeBoundingMappings();
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());
  filteredVertexPositions = filterMesh(seed, filteredMesh);
  seed.clear();
  seed.addMesh(filteredMesh);
  clearBoundingMappings();
}

void PreFilterPostFilterDecomposition:: feedback(
  mesh::Mesh& seed,
  std::map<int,std::vector<int> >& boundingVertexDistribution,
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

    //we need to merge the 2 filtering steps, each slave only holds local IDs
    std::vector<int> globalVertexIDs(numberOfVertices, -1);
    for (int i=0; i < numberOfVertices; i++) {
      globalVertexIDs[i] = boundingVertexDistribution[0][filteredVertexPositions[i]];
    }

    seed.getVertexDistribution()[0] = globalVertexIDs;

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      int numberOfVertices = -1;
      utils::MasterSlave::_communication->receive(numberOfVertices,rankSlave);
      std::vector<int> slaveVertexIDs(numberOfVertices,-1);
      if (numberOfVertices!=0) {
        utils::MasterSlave::_communication->receive(slaveVertexIDs.data(),numberOfVertices,rankSlave);
      }

      //we need to merge the 2 filtering steps, each slave only holds local IDs
      std::vector<int> globalVertexIDs(numberOfVertices,-1);
      for(int i=0;i<numberOfVertices;i++){
        globalVertexIDs[i] = boundingVertexDistribution[rankSlave][slaveVertexIDs[i]];
      }
      seed.getVertexDistribution()[rankSlave] = globalVertexIDs;
    }
  }

}



}}} // namespace precice, geometry, impl
