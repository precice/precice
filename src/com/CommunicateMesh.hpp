#pragma once

#include <string>
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace mesh {
class Mesh;
} // namespace mesh

namespace com {

/// Copies a Mesh object from a sender to a receiver.
class CommunicateMesh {
public:
  /// Constructor, takes communication to be used in transfer.
  explicit CommunicateMesh(
      com::PtrCommunication communication);

  /// Sends a constructed mesh to the receiver with given rank.
  void sendMesh(
      const mesh::Mesh &mesh,
      int               rankReceiver);

  /// Receives a mesh from the sender with given rank. Adds received mesh to mesh.
  void receiveMesh(
      mesh::Mesh &mesh,
      int         rankSender);

  void broadcastSendMesh(
      const mesh::Mesh &mesh);

  void broadcastReceiveMesh(
      mesh::Mesh &mesh);

private:
  logging::Logger _log{"com::CommunicateMesh"};

  /// Communication means used for the transfer of the geometry.
  com::PtrCommunication _communication;
};
} // namespace com
} // namespace precice
