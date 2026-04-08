#include "com/Extra.hpp"

#include "com/SerializedMesh.hpp"
#include "com/SerializedPartitioning.hpp"

namespace precice::com {

void sendMesh(IntraCommunication &communication, int rankReceiver, const mesh::Mesh &mesh)
{
  serialize::SerializedMesh::serialize(mesh).send(communication, rankReceiver);
}

void receiveMesh(IntraCommunication &communication, int rankSender, mesh::Mesh &mesh)
{
  serialize::SerializedMesh::receive(communication, rankSender).addToMesh(mesh);
}

void broadcastSendMesh(IntraCommunication &communication, const mesh::Mesh &mesh)
{
  serialize::SerializedMesh::serialize(mesh).broadcastSend(communication);
}

void broadcastReceiveMesh(IntraCommunication &communication, mesh::Mesh &mesh)
{
  serialize::SerializedMesh::broadcastReceive(communication).addToMesh(mesh);
}

void sendConnectionMap(IntraCommunication &communication, int rankReceiver, const mesh::Mesh::ConnectionMap &cm)
{
  serialize::SerializedConnectionMap::serialize(cm).send(communication, rankReceiver);
}

void receiveConnectionMap(IntraCommunication &communication, int rankSender, mesh::Mesh::ConnectionMap &cm)
{
  cm = serialize::SerializedConnectionMap::receive(communication, rankSender).toConnectionMap();
}

void broadcastSendConnectionMap(IntraCommunication &communication, const mesh::Mesh::ConnectionMap &cm)
{
  serialize::SerializedConnectionMap::serialize(cm).broadcastSend(communication);
}

void broadcastReceiveConnectionMap(IntraCommunication &communication, mesh::Mesh::ConnectionMap &cm)
{
  cm = serialize::SerializedConnectionMap::broadcastReceive(communication).toConnectionMap();
}

void sendBoundingBox(IntraCommunication &communication, int rankReceiver, const mesh::BoundingBox &bb)
{
  serialize::SerializedBoundingBox::serialize(bb).send(communication, rankReceiver);
}

void receiveBoundingBox(IntraCommunication &communication, int rankSender, mesh::BoundingBox &bb)
{
  bb = serialize::SerializedBoundingBox::receive(communication, rankSender).toBoundingBox();
}

void sendBoundingBoxMap(IntraCommunication &communication, int rankReceiver, const mesh::Mesh::BoundingBoxMap &bbm)
{
  serialize::SerializedBoundingBoxMap::serialize(bbm).send(communication, rankReceiver);
}

void receiveBoundingBoxMap(IntraCommunication &communication, int rankSender, mesh::Mesh::BoundingBoxMap &bbm)
{
  bbm = serialize::SerializedBoundingBoxMap::receive(communication, rankSender).toBoundingBoxMap();
}

void broadcastSendBoundingBoxMap(IntraCommunication &communication, const mesh::Mesh::BoundingBoxMap &bbm)
{
  serialize::SerializedBoundingBoxMap::serialize(bbm).broadcastSend(communication);
}

void broadcastReceiveBoundingBoxMap(IntraCommunication &communication, mesh::Mesh::BoundingBoxMap &bbm)
{
  bbm = serialize::SerializedBoundingBoxMap::broadcastReceive(communication).toBoundingBoxMap();
}

} // namespace precice::com
