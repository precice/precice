#include "partition/ProvidedBoundingBox.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EventTimings.hpp"
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
  TRACE();

  if (_hasToSend) {

    // each rank sends its bb to master
    if (utils::MasterSlave::_slaveMode) { //slave
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).sendBoundingBox(_mesh->getBoundingBox(), 0);
    } else { // Master

      assertion(utils::MasterSlave::_rank == 0);
      assertion(utils::MasterSlave::_size > 1);

      // to store the collection of bounding boxes
      mesh::Mesh::BoundingBoxMap bbm;

      // initialize bbm with dummy data
      mesh::Mesh::BoundingBox initialBB;
      for (int i = 0; i < _dimensions; i++) {
        initialBB.push_back(std::make_pair(-1, -1));
      }
      for (int rank = 0; rank < utils::MasterSlave::_size; rank++) {
        bbm[rank] = initialBB;
      }

      // master stores its bb into bbm
      bbm[0] = _mesh->getBoundingBox();

      // master receives bbs from slaves and stores them in bbm
      for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
        com::CommunicateBoundingBox(utils::MasterSlave::_communication).receiveBoundingBox(bbm[rankSlave], rankSlave);
      }

      // master sends number of ranks and bbm to the other master
      _m2n->getMasterCommunication()->send(utils::MasterSlave::_size, 0);
      com::CommunicateBoundingBox(_m2n->getMasterCommunication()).sendBoundingBoxMap(bbm, 0);
    }
  }
}

void ProvidedBoundingBox::computeBoundingBox()
{
  TRACE();

  // size of the feedbackmap that is received here
  int remoteConnectionMapSize = 0;

  std::map<int, std::vector<int>> remoteConnectionMap;

  if (not utils::MasterSlave::_slaveMode) { //Master
    assertion(utils::MasterSlave::_size > 1);

    // master receives feedback map (map of other participant ranks -> connected ranks at this participant)
    // from other participants master
    _m2n->getMasterCommunication()->receive(remoteConnectionMapSize, 0);
    for (int i = 0; i < remoteConnectionMapSize; i++) {
      remoteConnectionMap[i] = {-1};
    }
    if (remoteConnectionMapSize != 0)
      com::CommunicateBoundingBox(_m2n->getMasterCommunication()).receiveConnectionMap(remoteConnectionMap, 0);

    // broadcast the received feedbackMap
    utils::MasterSlave::_communication->broadcast(remoteConnectionMapSize);
    if (remoteConnectionMapSize != 0) {
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendConnectionMap(remoteConnectionMap);
    }

    // master checks which ranks are connected to it
    for (auto &remoteRank : remoteConnectionMap) {
      for (auto &includedRank : remoteRank.second) {
        if (utils::MasterSlave::_rank == includedRank) {
          _mesh->getConnectedRanks().push_back(remoteRank.first);
        }
      }
    }

  } else { // Slave

    utils::MasterSlave::_communication->broadcast(remoteConnectionMapSize, 0);

    for (int i = 0; i < remoteConnectionMapSize; i++) {
      remoteConnectionMap[i] = {-1};
    }
    if (remoteConnectionMapSize != 0)
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveConnectionMap(remoteConnectionMap);

    for (auto &remoteRank : remoteConnectionMap) {
      for (auto &includedRanks : remoteRank.second) {
        if (utils::MasterSlave::_rank == includedRanks) {
          _mesh->getConnectedRanks().push_back(remoteRank.first);
        }
      }
    }
  }
}

// these functions will be implemented in package 3
void ProvidedBoundingBox::communicate()
{
}

void ProvidedBoundingBox::compute()
{
}

void ProvidedBoundingBox::createOwnerInformation()
{
}

} // namespace partition
} // namespace precice
