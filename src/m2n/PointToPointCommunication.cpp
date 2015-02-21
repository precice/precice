// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "PointToPointCommunication.hpp"

#include "com/Request.hpp"
#include "com/SharedPointer.hpp"
#include "utils/MasterSlave.hpp"
#include "mesh/Mesh.hpp"

#include <vector>

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

// The complexity of `scatter' function is (total number of indices in
// `thisVertexDistribution') * (total number of indices in
// `otherVertexDistribution').
void
scatter(com::PtrCommunication communication,
        // `m' is output communication map for the master process.
        std::map<int, std::vector<int>>& m,
        // `thisVertexDistribution' is input vertex distribution from this
        // participant.
        std::map<int, std::vector<int>> const& thisVertexDistribution,
        // `otherVertexDistribution' is input vertex distribution from other
        // participant.
        std::map<int, std::vector<int>> const& otherVertexDistribution,
        std::function<void(int, std::map<int, std::vector<int>> const&)>
            function = [](int r, std::map<int, std::vector<int>> const& m) {}) {
  std::map<int, std::vector<int>> communicationMap;
  int i;

  forMapOfRanges(thisVertexDistribution,
                 [&](int thisRank) mutable { i = 0; },
                 [&](int thisRank, int thisIndex) mutable {
                   forMapOfRanges(otherVertexDistribution,
                                  [=, &communicationMap](
                                      int otherRank, int otherIndex) mutable {
                     if (thisIndex == otherIndex)
                       communicationMap[otherRank].push_back(i);
                   });

                   i++;
                 },
                 [&](int thisRank) mutable {
                   function(thisRank, communicationMap);

                   if (thisRank == utils::MasterSlave::_rank)
                     m = std::move(communicationMap);
                   else
                     send(communication, communicationMap, thisRank);

                   communicationMap.clear();
                 });

  if (thisVertexDistribution.find(utils::MasterSlave::_rank) ==
      thisVertexDistribution.end())
    m.clear();

  for (int rank = 0; rank < utils::MasterSlave::_size; ++rank) {
    if (rank == utils::MasterSlave::_rank)
      continue;

    if (thisVertexDistribution.find(rank) == thisVertexDistribution.end())
      send(communication, communicationMap, rank);
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
PointToPointCommunication::acceptConnection(std::string const& nameAcceptor,
                                            std::string const& nameRequester) {
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);

  assertion(not _isConnected);

  // The current participant is the *acceptor* participant and its processes are
  // *acceptor* processes.
  _isAcceptor = true;

  if (utils::MasterSlave::_masterMode) {
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();

    c->acceptConnection(nameAcceptor, nameRequester, 0, 1);
    // -------------------------------------------------------------------------
    // Exchange ranks of participants' master processes.
    int requesterMasterRank;

    c->send(utils::MasterSlave::_masterRank, 0);
    c->receive(requesterMasterRank, 0);
    // -------------------------------------------------------------------------
    // Exchange vertex distributions.
    auto& vertexDistribution = _mesh->getVertexDistribution();
    std::map<int, std::vector<int>> requesterVertexDistribution;

    m2n::send(c, vertexDistribution, 0);
    m2n::receive(c, requesterVertexDistribution, 0);
    // -------------------------------------------------------------------------
    // Vector of point-to-point communicator sizes. Provies one-to-one mapping
    // between acceptor process ranks (in the current participant) and
    // point-to-point communicator sizes (to be used during point-to-point
    // connection establishment).
    std::vector<int> communicatorSizes(utils::MasterSlave::_size, 0);

    // Iteratively construct different instances of communication map (from the
    // two vertex distributions), each of which corresponds to a single acceptor
    // process (in the current participant), and scatter them across their
    // corresponding acceptor processes (including the master acceptor process).
    m2n::scatter(
        utils::MasterSlave::_communication,
        _communicationMap,
        vertexDistribution,
        requesterVertexDistribution,
        // NOTE:
        // This lambda callback is invoked for each communication map (after its
        // construction is finished). The `scatter' routine is implemented with
        // memory efficiency in mind, and, therefore, as soon as this callback
        // returns, the memory occupied by (the current) `communicationMap' is
        // freed (to accomodate the next one). Thus, in order to prevent
        // undefined behavior, keeping any references to (the current)
        // `communicationMap' after this callback returns is strongly
        // discouraged.
        [&communicatorSizes](
            int rank, std::map<int, std::vector<int>> const& communicationMap) {
          // Initialize `communicatorSizes'. Notice how each acceptor process
          // rank `rank' (in the current participant) corresponds to the size of
          // `communicationMap' for this acceptor process.
          communicatorSizes[rank] = communicationMap.size();
        });

    // Send `communicatorSizes' to the requester participant.
    m2n::send(c, communicatorSizes, 0);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    // Receive the corresponding (to the current acceptor process) communication
    // map that is being scattered by the master acceptor process.
    m2n::receive(utils::MasterSlave::_communication, _communicationMap, 0);
  }

// NOTE:
// Change 0 to 1 to print `_communicationMap'.
#if 0
  print(_communicationMap);
#endif

  if (_communicationMap.size() == 0) {
    _isConnected = true;

    return;
  }

  // Accept point-to-point connections between the current acceptor process (in
  // the current participant) with rank `utils::MasterSlave::_rank' and
  // (multiple) requester proccesses (in the requester participant).
  auto c = _communicationFactory->newCommunication();

  c->acceptConnection(nameAcceptor + std::to_string(utils::MasterSlave::_rank),
                      nameRequester,
                      0,
                      1);

  // NOTE:
  // On the acceptor participant side, the communication object `c' behaves as a
  // server, i.e. it implicitly accepts multiple connections to requester
  // processes (in the requester participant) according to point-to-point
  // communicator size (in `communicatorSizes') used during point-to-point
  // connection requests made by requester processes (in the requester
  // participant). As a result, only one communication object `c' is needed to
  // satisfy `_communicationMap', and, therefore, for data structure consistency
  // of `_communications' with the requester participant side, we simply
  // duplicate references to the same communication object `c'.
  for (auto const& i : _communicationMap) {
    int requesterRank = i.first;

    _communications[requesterRank] = c;
  }

  _isConnected = true;
}

void
PointToPointCommunication::requestConnection(std::string const& nameAcceptor,
                                             std::string const& nameRequester) {
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);

  assertion(not _isConnected);

  // The current participant is the *requester* participant and its processes
  // are *requester* processes.
  _isAcceptor = false;

  std::vector<int> communicationRanks;
  std::vector<int> communicatorSizes;

  if (utils::MasterSlave::_masterMode) {
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();

    c->requestConnection(nameAcceptor, nameRequester, 0, 1);
    // -------------------------------------------------------------------------
    // Exchange ranks of participants' master processes.
    int acceptorMasterRank;

    c->receive(acceptorMasterRank, 0);
    c->send(utils::MasterSlave::_masterRank, 0);
    // -------------------------------------------------------------------------
    // Exchange vertex distributions.
    auto& vertexDistribution = _mesh->getVertexDistribution();
    std::map<int, std::vector<int>> acceptorVertexDistribution;

    m2n::receive(c, acceptorVertexDistribution, 0);
    m2n::send(c, vertexDistribution, 0);
    // -------------------------------------------------------------------------
    // Iteratively construct different instances of communication map (from the
    // two vertex distributions), each of which corresponds to a single
    // requester process (in the current participant), and scatter them across
    // their corresponding requester processes (including the master requester
    // process).
    m2n::scatter(utils::MasterSlave::_communication,
                 _communicationMap,
                 vertexDistribution,
                 acceptorVertexDistribution);

    // Receive `communicatorSizes' from the acceptor participant.
    m2n::receive(c, communicatorSizes, 0);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    // Receive the corresponding (to the current requester process)
    // communication map that is being scattered by the master requester
    // process.
    m2n::receive(utils::MasterSlave::_communication, _communicationMap, 0);
  }

// NOTE:
// Change 0 to 1 to print `_communicationMap'.
#if 0
  print(_communicationMap);
#endif

  //   ┌───┐  Legend
  //   ↓   │ ┌──────────────────────────────────────────────────────────────┐
  // ┌─┴─┐ │ │ ↓ — send/receive of `communicationRanks' and                 │
  // │ 0 │ │ │     `communicatorSizes';                                     │
  // └─┬─┘ │ │ r — global process rank in the requester participant;        │
  //   ↓   │ │ s — global communicator size in the requester participant.   │
  // ┌─┴─┐ │ └──────────────────────────────────────────────────────────────┘
  // │ 1 │ │
  // └─┬─┘ │
  //   ↓   │
  //   ⁞   │
  //   ↓   │
  // ┌─┴─┐ │
  // │ r │ │
  // └─┬─┘ │
  //   ↓   │
  //   ⁞   │
  //   ↓   │
  // ┌─┴─┐ │
  // │s-1│ │
  // └─┬─┘ │
  //   ↓   │
  //   └───┘
  //
  // Figure 1. Point-to-point connection request scheme

  // Point-to-point connection request scheme (Figure 1) is strictly sequential
  // in a sense that only one requester process (in the current participant) at
  // time requests point-to-point connection to (multiple) acceptor proccesses
  // (in the acceptor participant). This is due to how corresponding
  // point-to-point communication ranks (in `communicationRanks') have to be
  // incremented and sent further (to the next requester process). The execution
  // path starts and ends in the master requester process (in the current
  // participant).

  if (utils::MasterSlave::_masterMode) {
    // As discussed previously, the execution path of the point-to-point
    // connection request scheme (Figure 1) starts in the master requester
    // process (here). Thus, initialize `communicationRanks' and proceed.
    communicationRanks.resize(communicatorSizes.size(), -1);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    // All other processes should block until they receive (one by one) properly
    // incremented versions of `communicationRanks'.
    receivePrevious(_communicationFactory, communicationRanks);
    receivePrevious(_communicationFactory, communicatorSizes);

    assertion(communicationRanks.size() == communicatorSizes.size());
  }

  // Request point-to-point connections between the current requester process
  // (in the current participant) and (multiple) acceptor proccesses (in the
  // acceptor participant) with ranks `acceptorRank' according to communication
  // map.
  for (auto const& i : _communicationMap) {
    int acceptorRank = i.first;

    // Pay attention to how point-to-point communication rank (in
    // `communicationRanks') corresponding to `acceptorRank' is pre-incremented.
    int rank = ++communicationRanks[acceptorRank];
    int size = communicatorSizes[acceptorRank];

    _communications[acceptorRank] = _communicationFactory->newCommunication();

    _communications[acceptorRank]->requestConnection(
        nameAcceptor + std::to_string(acceptorRank), nameRequester, rank, size);

    // NOTE:
    // On the requester participant side, the communication objects behave as
    // clients, i.e. each of them requests only one connection to acceptor
    // process (in the acceptor participant).
  }

  // Send new incremented version of `communicationRanks' to the next requester
  // process (Figure 1).
  sendNext(_communicationFactory, communicationRanks);
  sendNext(_communicationFactory, communicatorSizes);

  if (utils::MasterSlave::_masterMode) {
    // As discussed previously, the execution path of the point-to-point
    // connection request scheme (Figure 1) ends in the master requester process
    // (here). Thus, receive the last incremented version of
    // `communicationRanks'.
    receivePrevious(_communicationFactory, communicationRanks);
    receivePrevious(_communicationFactory, communicatorSizes);

    assertion(communicationRanks.size() == communicatorSizes.size());

    // Perform an important sanity check for the last incremented version of
    // `communicationRanks'.
    for (int acceptorRank = 0;
         acceptorRank <
             std::min(communicationRanks.size(), communicatorSizes.size());
         ++acceptorRank) {
      preciceCheck(communicationRanks[acceptorRank] ==
                       communicatorSizes[acceptorRank] - 1,
                   "requestConnection()",
                   "Local communication rank/size inconsistency!");
    }
  }

  _isConnected = true;
}

void
PointToPointCommunication::closeConnection() {
  // TODO
  _isConnected = false;
}

void
PointToPointCommunication::send(double* itemsToSend,
                                int size,
                                int valueDimension) {
  assertion(_communicationMap.size() == _communications.size());

  if (_communicationMap.size() == 0) {
    preciceCheck(
        size == 0, "send()", "Can't send anything from disconnected process!");

    return;
  }

  std::vector<com::PtrRequest> requests;

  requests.reserve(_communicationMap.size());

  int rank = 0;

  for (auto i : _communicationMap) {
    auto receiverRank = i.first;
    auto indices = i.second;

    auto c = _communications[receiverRank];

    std::vector<double> data;

    data.reserve(indices.size() * valueDimension);

    for (auto index : indices) {
      for (int d = 0; d < valueDimension; ++d) {
        data.push_back(itemsToSend[index * valueDimension + d]);
      }
    }

    auto request = c->aSend(data.data(), data.size(), rank);

    requests.push_back(request);

    if (_isAcceptor)
      rank++;
  }

  for (auto request : requests) {
    request->wait();
  }
}

void
PointToPointCommunication::receive(double* itemsToReceive,
                                   int size,
                                   int valueDimension) {
  assertion(_communicationMap.size() == _communications.size());

  if (_communicationMap.size() == 0) {
    preciceCheck(size == 0,
                 "receive()",
                 "Can't receive anything to disconnected process!");

    return;
  }

  std::fill(itemsToReceive, itemsToReceive + size, 0);

  int rank = 0;

  for (auto i : _communicationMap) {
    auto otherRank = i.first;
    auto indices = i.second;

    auto c = _communications[otherRank];

    std::vector<double> data;

    data.resize(indices.size() * valueDimension);

    c->receive(data.data(), data.size(), rank);

    {
      int i = 0;

      for (auto index : indices) {
        for (int d = 0; d < valueDimension; ++d) {
          itemsToReceive[index * valueDimension + d] +=
              data[i * valueDimension + d];
        }

        i++;
      }
    }

    if (_isAcceptor)
      rank++;
  }
}
}
} // namespace precice, m2n
