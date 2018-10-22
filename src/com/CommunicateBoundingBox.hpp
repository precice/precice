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
      int                            rankReceiver);

  void receiveBoundingBox(
      mesh::Mesh::BoundingBox &bb,
      int                      rankSender);

  void sendBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &bbm,
      int                         rankReceiver);

  void receiveBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &bbm,
      int                         rankSender);

  void sendFeedbackMap(
      mesh::Mesh::FeedbackMap &fbm,
      int                         rankReceiver);

  void receiveFeedbackMap(
      mesh::Mesh::FeedbackMap &fbm,
      int                         rankSender);

  /// This method sends the set of bounding boxes (gathered in the master rank) to the other particpants master.
  /// Here we assume that the receiving particpant already knows about the size of sending participant's communicator.
  void broadcastSendBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &bbm);

  /// Receiveing participant's master calls this method to receive the set of bounding boxes sent by other
  /// participant's master. 
  void broadcastReceiveBoundingBoxMap(
      mesh::Mesh::BoundingBoxMap &bbm);

  void broadcastSendFeedbackMap(
    mesh::Mesh::FeedbackMap &fbm);

  void broadcastReceiveFeedbackMap(
    mesh::Mesh::FeedbackMap &fbm);
  
private:
  logging::Logger _log{"com::CommunicateBoundingBox"};

  /// Communication means used for the transfer of the geometry.
  com::PtrCommunication _communication;
};
} // namespace com
} // namespace precice
