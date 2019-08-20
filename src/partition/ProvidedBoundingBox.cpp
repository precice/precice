#include "partition/ProvidedBoundingBox.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EventUtils.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Parallel.hpp"

using precice::utils::Event;

namespace precice
{
namespace partition
{

ProvidedBoundingBox::ProvidedBoundingBox(mesh::PtrMesh mesh,
                                         bool          hasToSend,
                                         double        safetyFactor)
    : Partition(mesh),
      _hasToSend(hasToSend),
      _dimensions(mesh->getDimensions()),
      _safetyFactor(safetyFactor)
{
}

void ProvidedBoundingBox::communicateBoundingBox()
{
  PRECICE_TRACE();

  if (!_hasToSend)
    return;

  // each rank sends its bb to master
  if (utils::MasterSlave::isSlave()) { //slave
    com::CommunicateBoundingBox(utils::MasterSlave::_communication).sendBoundingBox(_mesh->getBoundingBox(), 0);
  } else { // Master

    PRECICE_ASSERT(utils::MasterSlave::getRank() == 0);
    PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

    // to store the collection of bounding boxes
    mesh::Mesh::BoundingBoxMap bbm;

    // initialize bbm with dummy data
    mesh::Mesh::BoundingBox initialBB;
    for (int i = 0; i < _dimensions; i++) {
      initialBB.push_back(std::make_pair(-1, -1));
    }
    for (int rank = 0; rank < utils::MasterSlave::getSize(); rank++) {
      bbm[rank] = initialBB;
    }

    // master stores its bb into bbm
    bbm[0] = _mesh->getBoundingBox();

    // master receives bbs from slaves and stores them in bbm
    for (int rankSlave = 1; rankSlave < utils::MasterSlave::getSize(); rankSlave++) {
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).receiveBoundingBox(bbm[rankSlave], rankSlave);
    }

    // master sends number of ranks and bbm to the other master
    _m2ns[0]->getMasterCommunication()->send(utils::MasterSlave::getSize(), 0);
    com::CommunicateBoundingBox(_m2ns[0]->getMasterCommunication()).sendBoundingBoxMap(bbm, 0);
  }
}

void ProvidedBoundingBox::computeBoundingBox()
{
  if (!_hasToSend)
    return;
  
  PRECICE_TRACE();

  // size of the feedbackmap
  int remoteConnectionMapSize = 0;
  std::vector<int> connectedRanksList;

  std::map<int, std::vector<int>> remoteConnectionMap;

  if (not utils::MasterSlave::isSlave()) { //Master
    PRECICE_ASSERT(utils::MasterSlave::getSize() > 1);

    // master receives feedback map (map of other participant ranks -> connected ranks at this participant)
    // from other participants master
    _m2ns[0]->getMasterCommunication()->receive(connectedRanksList, 0);
    remoteConnectionMapSize = connectedRanksList.size();
    
    for (auto &rank : connectedRanksList) {
      remoteConnectionMap[rank] = {-1};
    }
    if (remoteConnectionMapSize != 0)
      com::CommunicateBoundingBox(_m2ns[0]->getMasterCommunication()).receiveConnectionMap(remoteConnectionMap, 0);

    // broadcast the received feedbackMap
    utils::MasterSlave::_communication->broadcast(connectedRanksList);
    if (remoteConnectionMapSize != 0) {
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendConnectionMap(remoteConnectionMap);
    }

    // master checks which ranks are connected to it
    for (auto &remoteRank : remoteConnectionMap) {
      for (auto &includedRank : remoteRank.second) {
        if (utils::MasterSlave::getRank() == includedRank) {
          _mesh->getConnectedRanks().push_back(remoteRank.first);
        }
      }
    }

  } else { // Slave

    utils::MasterSlave::_communication->broadcast(connectedRanksList, 0);

    if (connectedRanksList.size() != 0)
    {
      for (auto &rank : connectedRanksList) {
        remoteConnectionMap[rank] = {-1};
      }
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveConnectionMap(remoteConnectionMap);
    }

    for (auto &remoteRank : remoteConnectionMap) {
      for (auto &includedRanks : remoteRank.second) {
        if (utils::MasterSlave::getRank() == includedRanks) {
          _mesh->getConnectedRanks().push_back(remoteRank.first);
        }
      }
    }
  }
}


void ProvidedBoundingBox::communicate()
{
  // each rank send its mesh partition to connected ranks at other participant
  _m2n->broadcastSendLocalMesh(*_mesh);  
}

void ProvidedBoundingBox::compute()
{}

void ProvidedBoundingBox::createOwnerInformation()
{}

} // namespace partition
} // namespace precice
