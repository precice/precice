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
  Decomposition ( dimensions ),
  _bb(),
  _safetyGap(0.0),
  _safetyFactor(safetyFactor),
  _filterByMapping(false)
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
  postFilter(seed, boundingVertexDistribution, filteredVertexPositions);
  feedback(seed, boundingVertexDistribution, filteredVertexPositions);
}

void PreFilterPostFilterDecomposition:: preFilter(
  mesh::Mesh& seed,
  std::map<int,std::vector<int> >& boundingVertexDistribution)
{
  preciceTrace1 ( "preFilter()", utils::MasterSlave::_rank );
  preciceInfo("preFilter()", "Pre-filter mesh " << seed.getName() );
  Event e("pre-filter");

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
    seed.setGlobalNumberOfVertices(seed.vertices().size());

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
  std::map<int,std::vector<int> >& boundingVertexDistribution,
  std::vector<int>& filteredVertexPositions)
{
  preciceTrace1 ( "postFilter()", utils::MasterSlave::_rank );
  preciceInfo("postFilter()", "Post-filter mesh " << seed.getName() );
  Event e("post-filter");

  seed.computeState();
  computeBoundingMappings();
  preciceInfo("scatterMesh()", "Filter mesh " << seed.getName() );
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, seed.isFlipNormals());
  filteredVertexPositions = filterMesh(seed, filteredMesh);
  seed.clear();
  seed.addMesh(filteredMesh);
  clearBoundingMappings();
}


void PreFilterPostFilterDecomposition:: mergeBoundingBoxes(mesh::Mesh::BoundingBox& bb){
  if (_boundingFromMapping.use_count()>0) {
    auto bb1 = _boundingFromMapping->getOutputMesh()->getBoundingBox();
    for (int d=0; d < _dimensions; d++) {
      if (bb[d].first > bb1[d].first) bb[d].first = bb1[d].first;
      if (bb[d].second < bb1[d].second) bb[d].second = bb1[d].second;
    }
  }
  if (_boundingToMapping.use_count()>0) {
    auto bb2 = _boundingToMapping->getInputMesh()->getBoundingBox();
    for (int d=0; d<_dimensions; d++) {
      if (bb[d].first > bb2[d].first) bb[d].first = bb2[d].first;
      if (bb[d].second < bb2[d].second) bb[d].second = bb2[d].second;
    }
  }
}

bool PreFilterPostFilterDecomposition:: doesVertexContribute(
  const mesh::Vertex& vertex)
{
  if (_filterByMapping) {
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
  else { //filter by bounding box
    for (int d=0; d<_dimensions; d++) {
      if (vertex.getCoords()[d] < _bb[d].first - _safetyGap || vertex.getCoords()[d] > _bb[d].second + _safetyGap) {
        return false;
      }
    }
    return true;
  }
}



}}} // namespace precice, geometry, impl
