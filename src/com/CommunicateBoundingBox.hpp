#pragma once
#include <map>
#include <string>
#include <vector>
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace mesh {
class BoundingBox;
} // namespace mesh

namespace com {

/// Copies either a bounding box around a mesh partition or complete maps of bounding boxes from a sender to a receiver.
class CommunicateBoundingBox {
public:
  /// Constructor, takes communication to be used in transfer.
  explicit CommunicateBoundingBox(
      com::PtrCommunication communication);

  void sendBoundingBox(
      const mesh::BoundingBox &bb,
      int                      rankReceiver);

  void receiveBoundingBox(
      mesh::BoundingBox &bb,
      int                rankSender);

  void sendBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &bbm,
      int                         rankReceiver);

  void receiveBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &bbm,
      int                         rankSender);

  void sendConnectionMap(
      std::map<int, std::vector<int>> const &fbm,
      int                                    rankReceiver);

  void receiveConnectionMap(
      std::map<int, std::vector<int>> &fbm,
      int                              rankSender);

  /// This method broadcasts the set of bounding boxes (gathered in the master rank) to the slaves.
  void broadcastSendBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &bbm);

  /// Slaves call this method to receive the set of bounding boxes sent by the master.
  void broadcastReceiveBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &bbm);

  void broadcastSendConnectionMap(
      std::map<int, std::vector<int>> const &fbm);

  void broadcastReceiveConnectionMap(
      std::map<int, std::vector<int>> &fbm);

private:
  logging::Logger _log{"com::CommunicateBoundingBox"};

  /// Communication means used for the transfer of the geometry.
  com::PtrCommunication _communication;
};
} // namespace com
} // namespace precice
