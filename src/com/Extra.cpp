#include "com/Extra.hpp"

#include "com/SerializedMesh.hpp"
#include "com/SerializedPartitioning.hpp"
#include "com/SerializedConnectionInfo.hpp"

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

void sendConnectionInfo(Communication &communication, int rankReceiver, const std::string &connectionInfo)
{
  communication.send(connectionInfo, rankReceiver);
}

void receiveConnectionInfo(Communication &communication, int rankSender, std::string &connectionInfo)
{
  communication.receive(connectionInfo, rankSender);
}

void sendConnectionInfoMap(Communication &communication, int rankReceiver, const serialize::SerializedConnectionInfoMap::ConnectionInfoMap &connectionInfoMap)
{
  serialize::SerializedConnectionInfoMap::serialize(connectionInfoMap).send(communication, rankReceiver);
}

void receiveConnectionInfoMap(Communication &communication, int rankSender, serialize::SerializedConnectionInfoMap::ConnectionInfoMap &connectionInfoMap)
{
  connectionInfoMap = serialize::SerializedConnectionInfoMap::receive(communication, rankSender).toConnectionInfoMap();
}

void broadcastSendConnectionInfoMap(Communication &communication, const serialize::SerializedConnectionInfoMap::ConnectionInfoMap &connectionInfoMap)
{
  serialize::SerializedConnectionInfoMap::serialize(connectionInfoMap).broadcastSend(communication);
}

void broadcastReceiveConnectionInfoMap(Communication &communication, serialize::SerializedConnectionInfoMap::ConnectionInfoMap &connectionInfoMap)
{
  connectionInfoMap = serialize::SerializedConnectionInfoMap::broadcastReceive(communication).toConnectionInfoMap();
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
