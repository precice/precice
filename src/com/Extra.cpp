#include "com/Extra.hpp"

#include "com/SerializedMesh.hpp"
#include "com/SerializedPartitioning.hpp"

namespace precice::com {

void sendMesh(Communication &communication, int rankReceiver, const mesh::Mesh &mesh)
{
  serialize::SerializedMesh::serialize(mesh).send(communication, rankReceiver);
}

void receiveMesh(Communication &communication, int rankSender, mesh::Mesh &mesh)
{
  serialize::SerializedMesh::receive(communication, rankSender).addToMesh(mesh);
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
  serialize::SerializedConnectionMap::serialize(cm).send(communication, rankReceiver);
}

void receiveConnectionMap(Communication &communication, int rankSender, mesh::Mesh::ConnectionMap &cm)
{
  cm = serialize::SerializedConnectionMap::receive(communication, rankSender).toConnectionMap();
}

void broadcastSendConnectionMap(Communication &communication, const mesh::Mesh::ConnectionMap &cm)
{
  serialize::SerializedConnectionMap::serialize(cm).broadcastSend(communication);
}

void broadcastReceiveConnectionMap(Communication &communication, mesh::Mesh::ConnectionMap &cm)
{
  cm = serialize::SerializedConnectionMap::broadcastReceive(communication).toConnectionMap();
}

void sendBoundingBox(Communication &communication, int rankReceiver, const mesh::BoundingBox &bb)
{
  serialize::SerializedBoundingBox::serialize(bb).send(communication, rankReceiver);
}

void receiveBoundingBox(Communication &communication, int rankSender, mesh::BoundingBox &bb)
{
  bb = serialize::SerializedBoundingBox::receive(communication, rankSender).toBoundingBox();
}

void sendBoundingBoxMap(Communication &communication, int rankReceiver, const mesh::Mesh::BoundingBoxMap &bbm)
{
  serialize::SerializedBoundingBoxMap::serialize(bbm).send(communication, rankReceiver);
}

void receiveBoundingBoxMap(Communication &communication, int rankSender, mesh::Mesh::BoundingBoxMap &bbm)
{
  bbm = serialize::SerializedBoundingBoxMap::receive(communication, rankSender).toBoundingBoxMap();
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
