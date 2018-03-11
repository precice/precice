#pragma once
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"

namespace precice
{
namespace com
{

/// Copies either a bounding box around a mesh partition or complete maps of bounding boxes from a sender to a receiver.
class CommunicateBoundingBox
{
public:
  /// Constructor, takes communication to be used in transfer.
  explicit CommunicateBoundingBox(
    com::PtrCommunication communication);

  void sendBoundingBox(
    const mesh::Mesh::BoundingBox &bb,
    int rankReceiver);

  void receiveBoundingBox(
    mesh::Mesh::BoundingBox &bb,
    int rankSender);

  void sendBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm,
    int rankReceiver);

  void receiveBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm,
    int rankSender);

  void broadcastSendBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm);

  void broadcastReceiveBoundingBoxMap(
    mesh::Mesh::BoundingBoxMap &bbm);
  
private:
  static logging::Logger _log;

  /// Communication means used for the transfer of the geometry.
  com::PtrCommunication _communication;
};
}
} // namespace precice, com
