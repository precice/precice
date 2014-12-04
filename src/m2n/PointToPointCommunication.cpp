// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "PointToPointCommunication.hpp"
#include "com/Communication.hpp"
#include "com/SocketCommunication.hpp"
#include "utils/MasterSlave.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace m2n {

template <class Range, class Function>
void
forRange(Range& range, Function function) {
  for (auto& i : range) {
    function(i);
  }
}

template <class Map, class Function>
void
forMapOfRanges(Map& map, Function function) {
  for (auto& i : map) {
    for (auto& j : i.second) {
      function(i.first, j);
    }
  }
}

template <class Map,
          class BeforeRangeFunction,
          class RangeFunction,
          class AfterRangeFunction>
void
forMapOfRanges(Map& map,
               BeforeRangeFunction beforeRangeFunction,
               RangeFunction rangeFunction,
               AfterRangeFunction afterRangeFunction) {
  for (auto& i : map) {
    beforeRangeFunction(i.first);
    for (auto& j : i.second) {
      rangeFunction(i.first, j);
    }
    afterRangeFunction(i.first);
  }
}

void
send(com::PtrCommunication communication,
     std::map<int, std::vector<int>> const& m,
     int rankReceiver) {
  communication->send(static_cast<int>(m.size()), rankReceiver);

  for (auto& i : m) {
    auto& rank = i.first;
    auto& indices = i.second;

    communication->send(rank, rankReceiver);
    communication->send(static_cast<int>(indices.size()), rankReceiver);
    communication->send(
        const_cast<int*>(&indices[0]), indices.size(), rankReceiver);
  }
}

void
receive(com::PtrCommunication communication,
        std::map<int, std::vector<int>>& m,
        int rankSender) {
  m.clear();

  int size = 0;

  communication->receive(size, rankSender);

  while (size--) {
    int rank = 0;
    int size = 0;

    communication->receive(rank, rankSender);
    communication->receive(size, rankSender);

    m[rank].resize(size);

    communication->receive(&m[rank][0], size, rankSender);
  }
}

void
scatter(com::PtrCommunication communication,
        std::map<int, std::vector<int>>& m,
        std::map<int, std::vector<int>> const& thisVertexDistribution,
        std::map<int, std::vector<int>> const& otherVertexDistribution) {
  std::map<int, std::vector<int>> senderMap;
  int i;

  forMapOfRanges(thisVertexDistribution,
                 [&](int thisRank) mutable { i = 0; },
                 [&](int thisRank, int thisIndex) mutable {
                   forMapOfRanges(
                       otherVertexDistribution,
                       [=, &senderMap](int otherRank, int otherIndex) mutable {
                         if (thisIndex == otherIndex)
                           senderMap[otherRank].push_back(i);
                       });

                   i++;
                 },
                 [&](int thisRank) mutable {
                   if (thisRank == utils::MasterSlave::_rank)
                     m = std::move(senderMap);
                   else
                     send(communication, senderMap, thisRank);

                   senderMap.clear();
                 });

  if (thisVertexDistribution.find(utils::MasterSlave::_rank) ==
      thisVertexDistribution.end())
    m.clear();

  for (int rank = 0; rank < utils::MasterSlave::_size; ++rank) {
    if (rank == utils::MasterSlave::_rank)
      continue;

    if (thisVertexDistribution.find(rank) == thisVertexDistribution.end())
      send(communication, senderMap, rank);
  }
}

void
print(std::map<int, std::vector<int>> const& m) {
  std::ostringstream oss;

  oss << "rank: " << utils::MasterSlave::_rank << "\n";

  forMapOfRanges(m, [&oss](int rank, int index) {
    oss << rank << ":"
        << " " << index << "\n";
  });

  if (utils::MasterSlave::_masterMode) {
    std::string s;

    for (int rank = 1; rank < utils::MasterSlave::_size; ++rank) {
      utils::MasterSlave::_communication->receive(s, rank);

      oss << s;
    }

    std::cout << oss.str();
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    utils::MasterSlave::_communication->send(oss.str(), 0);
  }
}

tarch::logging::Log PointToPointCommunication::_log(
    "precice::m2n::PointToPointCommunication");

PointToPointCommunication::PointToPointCommunication(mesh::PtrMesh pMesh)
    : _pMesh(pMesh), _isConnected(false) {
}

PointToPointCommunication::~PointToPointCommunication() {
  if (isConnected()) {
    closeConnection();
  }
}

bool
PointToPointCommunication::isConnected() {
  return _isConnected;
}

void
PointToPointCommunication::acceptConnection(const std::string& nameAcceptor,
                                            const std::string& nameRequester,
                                            int acceptorProcessRank,
                                            int acceptorCommunicatorSize) {
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);

  if (utils::MasterSlave::_masterMode) {
    _communications[0] =
        com::PtrCommunication(new com::SocketCommunication("lo", 50000, ""));

    _communications[0]->acceptConnection(nameAcceptor, nameRequester, 0, 1);
    // -------------------------------------------------------------------------
    auto& vertexDistribution = _pMesh->getVertexDistribution();

    send(_communications[0], vertexDistribution, 0);
    // -------------------------------------------------------------------------
    std::map<int, std::vector<int>> requesterVertexDistribution;

    receive(_communications[0], requesterVertexDistribution, 0);
    // -------------------------------------------------------------------------
    scatter(utils::MasterSlave::_communication,
            _senderMap,
            vertexDistribution,
            requesterVertexDistribution);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    receive(utils::MasterSlave::_communication, _senderMap, 0);
  }

  print(_senderMap);
}

void
PointToPointCommunication::requestConnection(const std::string& nameAcceptor,
                                             const std::string& nameRequester,
                                             int requesterProcessRank,
                                             int requesterCommunicatorSize) {
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);

  if (utils::MasterSlave::_masterMode) {
    _communications[0] =
        com::PtrCommunication(new com::SocketCommunication("lo", 50000, ""));

    _communications[0]->requestConnection(nameAcceptor, nameRequester, 0, 1);
    // -------------------------------------------------------------------------
    std::map<int, std::vector<int>> acceptorVertexDistribution;

    receive(_communications[0], acceptorVertexDistribution, 0);
    // -------------------------------------------------------------------------
    auto& vertexDistribution = _pMesh->getVertexDistribution();

    send(_communications[0], vertexDistribution, 0);
    // -------------------------------------------------------------------------
    scatter(utils::MasterSlave::_communication,
            _senderMap,
            vertexDistribution,
            acceptorVertexDistribution);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    receive(utils::MasterSlave::_communication, _senderMap, 0);
  }

  print(_senderMap);
}

void
PointToPointCommunication::closeConnection() {
}

com::PtrCommunication
PointToPointCommunication::getMasterCommunication() {
  return _communications[0];
}

void
PointToPointCommunication::startSendPackage(int rankReceiver) {
}

void
PointToPointCommunication::finishSendPackage() {
}

int
PointToPointCommunication::startReceivePackage(int rankSender) {
  return -1;
}

void
PointToPointCommunication::finishReceivePackage() {
}

void
PointToPointCommunication::sendMaster(const std::string& itemToSend,
                                      int rankReceiver) {
}

void
PointToPointCommunication::sendMaster(int* itemsToSend,
                                      int size,
                                      int rankReceiver) {
}

void
PointToPointCommunication::sendMaster(double* itemsToSend,
                                      int size,
                                      int rankReceiver) {
}

void
PointToPointCommunication::sendMaster(double itemToSend, int rankReceiver) {
}

void
PointToPointCommunication::sendMaster(int itemToSend, int rankReceiver) {
}

void
PointToPointCommunication::sendMaster(bool itemToSend, int rankReceiver) {
}

int
PointToPointCommunication::receiveMaster(std::string& itemToReceive,
                                         int rankSender) {
  return -1;
}

int
PointToPointCommunication::receiveMaster(int* itemsToReceive,
                                         int size,
                                         int rankSender) {
  return -1;
}

int
PointToPointCommunication::receiveMaster(double* itemsToReceive,
                                         int size,
                                         int rankSender) {
  return -1;
}

int
PointToPointCommunication::receiveMaster(double& itemToReceive,
                                         int rankSender) {
  return -1;
}

int
PointToPointCommunication::receiveMaster(int& itemToReceive, int rankSender) {
  return -1;
}

int
PointToPointCommunication::receiveMaster(bool& itemToReceive, int rankSender) {
  return -1;
}

void
PointToPointCommunication::sendAll(utils::DynVector* itemsToSend,
                                   int size,
                                   int rankReceiver,
                                   mesh::PtrMesh mesh,
                                   int valueDimension) {
}

void
PointToPointCommunication::receiveAll(utils::DynVector* itemsToReceive,
                                      int size,
                                      int rankSender,
                                      mesh::PtrMesh mesh,
                                      int valueDimension) {
}

void
PointToPointCommunication::receiveAll(bool& itemToReceive, int rankSender) {
}

void
PointToPointCommunication::receiveAll(double& itemToReceive, int rankSender) {
}

void
PointToPointCommunication::sendAll(bool itemToSend, int rankReceiver) {
}

void
PointToPointCommunication::sendAll(double itemToSend, int rankReceiver) {
}
}
} // namespace precice, m2n
