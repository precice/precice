#pragma once

#include "SerializedConnectionInfo.hpp"
#include "com/Communication.hpp"
#include "mesh/Mesh.hpp"

namespace precice::com {

void sendMesh(Communication &communication, int rankReceiver, const mesh::Mesh &mesh);

void receiveMesh(Communication &communication, int rankSender, mesh::Mesh &mesh);

void broadcastSendMesh(Communication &communication, const mesh::Mesh &mesh);

void broadcastReceiveMesh(Communication &communication, mesh::Mesh &mesh);

void sendConnectionMap(Communication &communication, int rankReceiver, const mesh::Mesh::ConnectionMap &cm);

void receiveConnectionMap(Communication &communication, int rankSender, mesh::Mesh::ConnectionMap &cm);

void broadcastSendConnectionMap(Communication &communication, const mesh::Mesh::ConnectionMap &cm);

void broadcastReceiveConnectionMap(Communication &communication, mesh::Mesh::ConnectionMap &cm);

void sendBoundingBox(Communication &communication, int rankReceiver, const mesh::BoundingBox &bb);

void receiveBoundingBox(Communication &communication, int rankSender, mesh::BoundingBox &bb);

void sendBoundingBoxMap(Communication &communication, int rankReceiver, const mesh::Mesh::BoundingBoxMap &bbm);

void receiveBoundingBoxMap(Communication &communication, int rankSender, mesh::Mesh::BoundingBoxMap &bbm);

void sendConnectionInfo(Communication &communication, int rankReceiver, const std::string &connectionInfo);

void receiveConnectionInfo(Communication &communication, int rankSender, std::string &connectionInfo);

void sendConnectionInfoMap(Communication &communication, int rankReceiver, const serialize::SerializedConnectionInfoMap::ConnectionInfoMap &connectionInfoMap);

void receiveConnectionInfoMap(Communication &communication, int rankSender, serialize::SerializedConnectionInfoMap::ConnectionInfoMap &connectionInfoMap);

void broadcastSendConnectionInfoMap(Communication &communication, const serialize::SerializedConnectionInfoMap::ConnectionInfoMap &connectionInfoMap);

void broadcastReceiveConnectionInfoMap(Communication &communication, serialize::SerializedConnectionInfoMap::ConnectionInfoMap &connectionInfoMap);

void broadcastSendBoundingBoxMap(Communication &communication, const mesh::Mesh::BoundingBoxMap &bbm);

void broadcastReceiveBoundingBoxMap(Communication &communication, mesh::Mesh::BoundingBoxMap &bbm);

} // namespace precice::com
