#pragma once

#include <vector>

#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"

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

private:
  logging::Logger _log{"com::CommunicateBoundingBox"};

  /// Communication means used for the transfer of the geometry.
  com::PtrCommunication _communication;
};
} // namespace com
} // namespace precice
