#include "com/Extra.hpp"

#include "com/SerializeMesh.hpp"
#include "com/SerializePartitioning.hpp"

namespace precice::com {

void sendMesh(Communication &communication, int rankReceiver, const mesh::Mesh &mesh)
{
  serialize::SerializedMesh::serialize(mesh).sendTo(communication, rankReceiver);
}

void receiveMesh(Communication &communication, int rankSender, mesh::Mesh &mesh)
{
  serialize::SerializedMesh::receiveFrom(communication, rankSender).addToMesh(mesh);
}

void broadcastSendMesh(Communication &communication, const mesh::Mesh &mesh)
{
  serialize::SerializedMesh::serialize(mesh).broadcastSend(communication);
}

void broadcastReceiveMesh(Communication &communication, mesh::Mesh &mesh)
{
  serialize::SerializedMesh::broadcastReceive(communication).addToMesh(mesh);
}

void sendConnectionMap(Communication &communication, int rankReceiver, const mesh::Mesh::ConnectionMap &cm)
{
  serialize::SerializedConnectionMap::serialize(cm).sendTo(communication, rankReceiver);
}

void receiveConnectionMap(Communication &communication, int rankSender, mesh::Mesh::ConnectionMap &cm)
{
  cm = serialize::SerializedConnectionMap::receiveFrom(communication, rankSender).toConnectionMap();
}

void broadcastSendConnectionMap(Communication &communication, const mesh::Mesh::ConnectionMap &cm)
{
  serialize::SerializedConnectionMap::serialize(cm).broadcastSend(communication);
}

void broadcastReceiveConnectionMap(Communication &communication, mesh::Mesh::ConnectionMap &cm)
{
  cm = serialize::SerializedConnectionMap::broadcastReceive(communication).toConnectionMap();
}

void sendBoundingBoxMap(Communication &communication, int rankReceiver, const mesh::Mesh::BoundingBoxMap &bbm)
{
  serialize::SerializedBoundingBoxMap::serialize(bbm).sendTo(communication, rankReceiver);
}

void receiveBoundingBoxMap(Communication &communication, int rankSender, mesh::Mesh::BoundingBoxMap &bbm)
{
  bbm = serialize::SerializedBoundingBoxMap::receiveFrom(communication, rankSender).toBoundingBoxMap();
}

void broadcastSendBoundingBoxMap(Communication &communication, const mesh::Mesh::BoundingBoxMap &bbm)
{
  serialize::SerializedBoundingBoxMap::serialize(bbm).broadcastSend(communication);
}

void broadcastReceiveBoundingBoxMap(Communication &communication, mesh::Mesh::BoundingBoxMap &bbm)
{
  bbm = serialize::SerializedBoundingBoxMap::broadcastReceive(communication).toBoundingBoxMap();
}

} // namespace precice::com
