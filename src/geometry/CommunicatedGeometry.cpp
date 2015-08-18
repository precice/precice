#include "CommunicatedGeometry.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"
#include "utils/EventTimings.hpp"

using precice::utils::Event;

namespace precice {
namespace geometry {

tarch::logging::Log CommunicatedGeometry:: _log ( "precice::geometry::CommunicatedGeometry" );

CommunicatedGeometry:: CommunicatedGeometry
(
  const utils::DynVector& offset,
  const std::string&      accessor,
  const std::string&      provider,
  impl::PtrDecomposition  decomposition)
  :
  Geometry ( offset ),
  _accessorName ( accessor ),
  _providerName ( provider ),
  _receivers (),
  _decomposition(decomposition)
{
  preciceTrace2 ( "CommunicatedGeometry()", accessor, provider );
}

void CommunicatedGeometry:: addReceiver
(
  const std::string&     receiver,
  m2n::M2N::SharedPointer m2n)
{
  preciceTrace1 ( "addReceiver()", receiver );
  assertion ( m2n.get() != nullptr );
  preciceCheck ( ! utils::contained(receiver, _receivers),
                 "addReceiver()", "Receiver \"" << receiver
                 << "\" has been added already to communicated geometry!" );
  preciceCheck ( receiver != _providerName, "addReceiver()",
                 "Receiver \"" << receiver << "\" cannot be the same as "
                 << "provider in communicated geometry!" );
  _receivers[receiver] = m2n;
}

void CommunicatedGeometry:: specializedCreate
(
  mesh::Mesh& seed )
{
  preciceTrace1 ( "specializedCreate()", seed.getName() );
  preciceCheck ( _receivers.size() > 0, "specializedCreate()",
                 "No receivers specified for communicated geometry to create "
                 << "mesh \"" << seed.getName() << "\"!" );
  if ( _accessorName == _providerName ) {
    sendMesh(seed);
  }
  else if ( utils::contained(_accessorName, _receivers) ) {
    receiveMesh(seed);
  }
  else {
    preciceError( "specializedCreate()", "Participant \"" << _accessorName
                  << "\" uses a communicated geometry to create mesh \""
                  << seed.getName()
                  << "\" but is neither provider nor receiver!" );
  }
}


void CommunicatedGeometry:: sendMesh(
  mesh::Mesh& seed)
{
  preciceTrace1 ( "sendMesh()", utils::MasterSlave::_rank );
  // Temporary globalMesh such that the master also keeps his local mesh (seed)
  mesh::Mesh globalMesh(seed.getName(), seed.getDimensions(), seed.isFlipNormals());

  if( not utils::MasterSlave::_slaveMode ){
    globalMesh.addMesh(seed); //add local master mesh to global mesh
  }

  // Gather Mesh
  preciceInfo("sendMesh()", "Gather mesh " << seed.getName() );
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode ) {
    Event e("gather mesh");
    seed.getVertexOffsets().resize(utils::MasterSlave::_size);
    if (utils::MasterSlave::_slaveMode) {
      com::CommunicateMesh(utils::MasterSlave::_communication).sendMesh( seed, 0 );
      utils::MasterSlave::_communication->broadcast(seed.getVertexOffsets().data(),utils::MasterSlave::_size,0);
    }
    else{ // Master
      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);

      int numberOfVertices = 0;
      // Vertices of master mesh part do already exist
      for (size_t i = 0; i < seed.vertices().size(); i++) {
        seed.getVertexDistribution()[0].push_back(numberOfVertices);
        numberOfVertices++;
      }
      seed.getVertexOffsets()[0] = numberOfVertices;

      for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
        int vertexCount1 = globalMesh.vertices().size();
        com::CommunicateMesh(utils::MasterSlave::_communication).receiveMesh ( globalMesh, rankSlave);
        int vertexCount2 = globalMesh.vertices().size();
        int vertexCountDiff = vertexCount2 - vertexCount1;
        preciceDebug("Received sub-mesh, from slave: " << rankSlave <<", vertexCount: " << vertexCountDiff);
        for (int i = 0; i < vertexCountDiff; i++) {
          seed.getVertexDistribution()[rankSlave].push_back(numberOfVertices);
          numberOfVertices++;
        }
        seed.getVertexOffsets()[rankSlave] = numberOfVertices;
      }
      utils::MasterSlave::_communication->broadcast(seed.getVertexOffsets().data(),utils::MasterSlave::_size);
    }
  }

  // Send (global) Mesh
  preciceInfo("sendMesh()", "Send global mesh " << seed.getName());
  Event e("send global mesh");
  if (not utils::MasterSlave::_slaveMode) {
    preciceCheck ( globalMesh.vertices().size() > 0,
                   "specializedCreate()", "Participant \"" << _accessorName
                   << "\" provides an invalid (possibly empty) mesh \""
                   << globalMesh.getName() << "\"!" );
    for (auto &pair : _receivers) {
      com::CommunicateMesh(pair.second->getMasterCommunication()).sendMesh ( globalMesh, 0 );
    }
  }
}

void CommunicatedGeometry:: receiveMesh(
  mesh::Mesh& seed)
{
  preciceInfo("receiveMesh()", "Receive global mesh " << seed.getName() );
  if (not utils::MasterSlave::_slaveMode) {
    Event e("receive global mesh");
    assertion ( seed.vertices().size() == 0 );
    assertion ( utils::contained(_accessorName, _receivers) );
    m2n::M2N::SharedPointer m2n ( _receivers[_accessorName] );
    com::CommunicateMesh(m2n->getMasterCommunication()).receiveMesh ( seed, 0 );
  }
  if (utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode){
    seed.getVertexOffsets().resize(utils::MasterSlave::_size);
    _decomposition->decompose(seed);
    broadcastOwnerInformation(seed);
  }
}

void CommunicatedGeometry:: broadcastOwnerInformation(
  mesh::Mesh& seed)
{
  preciceTrace ( "sendMesh()");
  //send global index back to slave (needed for petsc)
  if (utils::MasterSlave::_slaveMode) {
    int numberOfVertices = seed.vertices().size();
    if (numberOfVertices!=0) {
      std::vector<int> globalIndices(numberOfVertices, -1);
      utils::MasterSlave::_communication->receive(globalIndices.data(),numberOfVertices,0);
      preciceDebug("My global indices: " << globalIndices);
      seed.setGlobalIndices(globalIndices);
    }
  }
  else { // Master
    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      auto globalIndices = seed.getVertexDistribution()[rankSlave];
      int numberOfVertices = globalIndices.size();
      if (numberOfVertices!=0) {
        utils::MasterSlave::_communication->send(globalIndices.data(),numberOfVertices,rankSlave);
      }
    }
    seed.setGlobalIndices(seed.getVertexDistribution()[0]);
  }

  //send owner information to slave (needed for petsc)
  if (utils::MasterSlave::_slaveMode) {
    int numberOfVertices = seed.vertices().size();
    if (numberOfVertices!=0) {
      std::vector<int> ownerVec(numberOfVertices, -1);
      utils::MasterSlave::_communication->receive(ownerVec.data(),numberOfVertices,0);
      preciceDebug("My owner information: " << ownerVec);
      seed.setOwnerInformation(ownerVec);
    }
  }
  else { // Master
    std::vector<int> globalOwnerVec(seed.getVertexOffsets()[utils::MasterSlave::_size-1],0);
    for (int rankSlave = 0; rankSlave < utils::MasterSlave::_size; rankSlave++){
      auto globalIndices = seed.getVertexDistribution()[rankSlave];
      int numberOfVertices = globalIndices.size();
      std::vector<int> ownerVec(numberOfVertices,0);
      for(int i=0;i<numberOfVertices;i++){
        if(globalOwnerVec[globalIndices[i]] == 0){
          ownerVec[i] = 1;
          globalOwnerVec[globalIndices[i]] = 1;
        }
      }
      if (numberOfVertices!=0) {
        if(rankSlave==0){ //master own data
          seed.setOwnerInformation(ownerVec);
        }
        else{
          utils::MasterSlave::_communication->send(ownerVec.data(),numberOfVertices,rankSlave);
        }
      }
    }
#   ifdef Debug
      for(int i=0;i<seed.getVertexOffsets()[utils::MasterSlave::_size-1];i++){
        if(globalOwnerVec[i]==0){
          preciceWarning("scatterMesh()", "The Vertex with global index " << i << " of mesh: " << seed.getName()
              << " was completely filtered out, since it has no influence on any mapping.")
        }
      }
#   endif

  }
}

}} // namespace precice, geometry
