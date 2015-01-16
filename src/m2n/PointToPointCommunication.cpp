// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_OMP

#include "PointToPointCommunication.hpp"

#include "utils/MasterSlave.hpp"
#include "mesh/Mesh.hpp"

#include <omp.h>

namespace precice {
namespace m2n {

template <class Map, class Function>
void
forMap(Map& map, Function function) {
  for (auto& i : map) {
    function(i.first, i.second);
  }
}

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
     std::vector<int> const& v,
     int rankReceiver) {
  communication->send(static_cast<int>(v.size()), rankReceiver);
  communication->send(const_cast<int*>(&v[0]), v.size(), rankReceiver);
}

void
receive(com::PtrCommunication communication,
        std::vector<int>& v,
        int rankSender) {
  v.clear();

  int size = 0;

  communication->receive(size, rankSender);
  v.resize(size);
  communication->receive(&v[0], size, rankSender);
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
    send(communication, indices, rankReceiver);
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

    communication->receive(rank, rankSender);
    receive(communication, m[rank], rankSender);
  }
}

void
scatter(com::PtrCommunication communication,
        std::map<int, std::vector<int>>& m,
        std::map<int, std::vector<int>> const& thisVertexDistribution,
        std::map<int, std::vector<int>> const& otherVertexDistribution,
        std::function<void(int, std::map<int, std::vector<int>> const&)>
            function = [](int r, std::map<int, std::vector<int>> const& m) {}) {
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
                   function(thisRank, senderMap);

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
sendNext(com::PtrCommunicationFactory communicationFactory,
         std::vector<int> const& v) {
  int rank = utils::MasterSlave::_rank;
  int nextRank = (rank + 1) % utils::MasterSlave::_size;

  auto c = communicationFactory->newCommunication();

  // NOTE:
  // If more than two participants are connected in the same directory, then
  // this might cause collision. Perhaps utilizing requester name as a prefix
  // would solve the problem.
  c->requestConnection("Receiver" + std::to_string(nextRank),
                       "Sender" + std::to_string(rank),
                       0,
                       1);

  send(c, v, 0);
}

void
receivePrevious(com::PtrCommunicationFactory communicationFactory,
                std::vector<int>& v) {
  int rank = utils::MasterSlave::_rank;
  int previousRank =
      (utils::MasterSlave::_size + rank - 1) % utils::MasterSlave::_size;

  auto c = communicationFactory->newCommunication();

  // NOTE:
  // If more than two participants are connected in the same directory, then
  // this might cause collision. Perhaps utilizing requester name as a prefix
  // would solve the problem.
  c->acceptConnection("Receiver" + std::to_string(rank),
                      "Sender" + std::to_string(previousRank),
                      0,
                      1);

  receive(c, v, 0);
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

PointToPointCommunication::PointToPointCommunication(
    com::PtrCommunicationFactory communicationFactory, mesh::PtrMesh mesh)
    : DistributedCommunication(mesh)
    , _communicationFactory(communicationFactory)
    , _isConnected(false)
    , _isAcceptor(false) {
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
                                            const std::string& nameRequester) {
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);

  _isAcceptor = true;

  if (utils::MasterSlave::_masterMode) {
    auto c = _communicationFactory->newCommunication();

    c->acceptConnection(nameAcceptor, nameRequester, 0, 1);
    // -------------------------------------------------------------------------
    c->send(utils::MasterSlave::_masterRank, 0);
    // -------------------------------------------------------------------------
    int requesterMasterRank;

    c->receive(requesterMasterRank, 0);
    // -------------------------------------------------------------------------
    auto& vertexDistribution = _mesh->getVertexDistribution();

    m2n::send(c, vertexDistribution, 0);
    // -------------------------------------------------------------------------
    std::map<int, std::vector<int>> requesterVertexDistribution;

    m2n::receive(c, requesterVertexDistribution, 0);
    // -------------------------------------------------------------------------
    std::vector<int> sizes(utils::MasterSlave::_size, 0);

    scatter(
        utils::MasterSlave::_communication,
        _senderMap,
        vertexDistribution,
        requesterVertexDistribution,
        [&sizes](int rank, std::map<int, std::vector<int>> const& senderMap) {
          sizes[rank] = senderMap.size();
        });

    m2n::send(c, sizes, 0);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    m2n::receive(utils::MasterSlave::_communication, _senderMap, 0);
  }

  print(_senderMap);

  if (_senderMap.size() == 0)
    return;

  auto c = _communicationFactory->newCommunication();

  c->acceptConnection(nameAcceptor + std::to_string(utils::MasterSlave::_rank),
                      nameRequester,
                      0,
                      1);

  for (auto const& i : _senderMap) {
    int requesterRank = i.first;

    _communications[requesterRank] = c;
  }

  _isConnected = true;
}

void
PointToPointCommunication::requestConnection(const std::string& nameAcceptor,
                                             const std::string& nameRequester) {
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);

  _isAcceptor = false;

  std::vector<int> acceptorRanks;
  std::vector<int> acceptorSizes;

  if (utils::MasterSlave::_masterMode) {
    auto c = _communicationFactory->newCommunication();

    c->requestConnection(nameAcceptor, nameRequester, 0, 1);
    // -------------------------------------------------------------------------
    int acceptorMasterRank;

    c->receive(acceptorMasterRank, 0);
    // -------------------------------------------------------------------------
    c->send(utils::MasterSlave::_masterRank, 0);
    // -------------------------------------------------------------------------
    std::map<int, std::vector<int>> acceptorVertexDistribution;

    m2n::receive(c, acceptorVertexDistribution, 0);
    // -------------------------------------------------------------------------
    auto& vertexDistribution = _mesh->getVertexDistribution();

    m2n::send(c, vertexDistribution, 0);
    // -------------------------------------------------------------------------
    scatter(utils::MasterSlave::_communication,
            _senderMap,
            vertexDistribution,
            acceptorVertexDistribution);

    m2n::receive(c, acceptorSizes, 0);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    m2n::receive(utils::MasterSlave::_communication, _senderMap, 0);
  }

  print(_senderMap);

  if (utils::MasterSlave::_masterMode) {
    acceptorRanks.resize(acceptorSizes.size(), -1);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    receivePrevious(_communicationFactory, acceptorRanks);
    receivePrevious(_communicationFactory, acceptorSizes);

    assertion(acceptorRanks.size() == acceptorSizes.size());
  }

  for (auto const& i : _senderMap) {
    int acceptorRank = i.first;

    int rank = ++acceptorRanks[acceptorRank]; // Attention!
    int size = acceptorSizes[acceptorRank];

    _communications[acceptorRank] = _communicationFactory->newCommunication();

    _communications[acceptorRank]->requestConnection(
        nameAcceptor + std::to_string(acceptorRank), nameRequester, rank, size);
  }

  sendNext(_communicationFactory, acceptorRanks);
  sendNext(_communicationFactory, acceptorSizes);

  if (utils::MasterSlave::_masterMode) {
    receivePrevious(_communicationFactory, acceptorRanks);
    receivePrevious(_communicationFactory, acceptorSizes);

    assertion(acceptorRanks.size() == acceptorSizes.size());

    for (int acceptorRank = 0;
         acceptorRank < std::min(acceptorRanks.size(), acceptorSizes.size());
         ++acceptorRank) {
      preciceCheck(
          acceptorRanks[acceptorRank] == acceptorSizes[acceptorRank] - 1,
          "requestConnection()",
          "Local communication rank/size inconsistency!");
    }
  }

  _isConnected = true;
}

void
PointToPointCommunication::closeConnection() {
}

void
PointToPointCommunication::send(double* itemsToSend,
                                int size,
                                int valueDimension) {
  assertion(_senderMap.size() == _communications.size());

  if (_senderMap.size() == 0) {
    preciceCheck(size == 0,
                 "sendAll()",
                 "Can't send anything from disconnected process!");

    return;
  }

  omp_set_dynamic(0);

#pragma omp parallel num_threads(_senderMap.size())
  {
    // TODO: Might be reasonable to prepare an array of iterators in order to
    // improve performance.
    auto indices = std::next(_senderMap.begin(), omp_get_thread_num())->second;
    auto c = std::next(_communications.begin(), omp_get_thread_num())->second;

    int rank = _isAcceptor ? omp_get_thread_num() : 0;

    for (auto index : indices) {
      c->send(itemsToSend[index], rank);
    }
  }
}

void
PointToPointCommunication::receive(double* itemsToReceive,
                                   int size,
                                   int valueDimension) {
  assertion(_senderMap.size() == _communications.size());

  int rank = 0;

  for (auto i : _senderMap) {
    auto otherRank = i.first;
    auto indices = i.second;

    auto c = _communications[otherRank];

    for (auto index : indices) {
      double item;

      c->receive(item, rank);

      itemsToReceive[index] += item;
    }

    if (_isAcceptor)
      rank++;
  }
}
}
} // namespace precice, m2n

#endif // not PRECICE_NO_OMP
