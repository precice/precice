#pragma once

#include "com/IntraCommunication.hpp"
#include "mesh/Mesh.hpp"

namespace precice::com {

void sendMesh(IntraCommunication &communication, int rankReceiver, const mesh::Mesh &mesh);

void receiveMesh(IntraCommunication &communication, int rankSender, mesh::Mesh &mesh);

void broadcastSendMesh(IntraCommunication &communication, const mesh::Mesh &mesh);

void broadcastReceiveMesh(IntraCommunication &communication, mesh::Mesh &mesh);

void sendConnectionMap(IntraCommunication &communication, int rankReceiver, const mesh::Mesh::ConnectionMap &cm);

void receiveConnectionMap(IntraCommunication &communication, int rankSender, mesh::Mesh::ConnectionMap &cm);

void broadcastSendConnectionMap(IntraCommunication &communication, const mesh::Mesh::ConnectionMap &cm);

void broadcastReceiveConnectionMap(IntraCommunication &communication, mesh::Mesh::ConnectionMap &cm);

void sendBoundingBox(IntraCommunication &communication, int rankReceiver, const mesh::BoundingBox &bb);

void receiveBoundingBox(IntraCommunication &communication, int rankSender, mesh::BoundingBox &bb);

void sendBoundingBoxMap(IntraCommunication &communication, int rankReceiver, const mesh::Mesh::BoundingBoxMap &bbm);

void receiveBoundingBoxMap(IntraCommunication &communication, int rankSender, mesh::Mesh::BoundingBoxMap &bbm);

void broadcastSendBoundingBoxMap(IntraCommunication &communication, const mesh::Mesh::BoundingBoxMap &bbm);

void broadcastReceiveBoundingBoxMap(IntraCommunication &communication, mesh::Mesh::BoundingBoxMap &bbm);

} // namespace precice::com
