#pragma once
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"

namespace precice
{
namespace com
{

/// Copies a Mesh object from a sender to a receiver.
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
    int rankSender,
    int dim);

  void sendBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &gbb,
      int rankReceiver,
      int dim);

  void receiveBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &gbb,
      int rankSender,
      int dim);

  void sendBoundingBoxFeedBack(
    std::vector<int> &feedback,
    int rankReceiver);

  void receiveBoundingBoxFeedBack(
    std::map<int, std::vector<int>> &feedbackMap,
    int rankSender);

  void broadcastSendBoundingBox(
    mesh::Mesh::BoundingBoxMap &gBB,
    int dim);

  void broadcastReceiveBoundingBox(
    mesh::Mesh::BoundingBoxMap &gBB,
    int dim);
  
private:
  static logging::Logger _log;

  /// Communication means used for the transfer of the geometry.
  com::PtrCommunication _communication;
};
}
} // namespace precice, com
