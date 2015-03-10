// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "PointToPointCommunication.hpp"

#include "com/Request.hpp"
#include "com/SharedPointer.hpp"
#include "utils/EventTimings.hpp"
#include "utils/MasterSlave.hpp"
#include "mesh/Mesh.hpp"

#include <vector>

using precice::utils::Event;

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
        std::map<int, std::vector<int>> const& otherVertexDistribution) {
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
    , _isConnected(false) {
}

PointToPointCommunication::~PointToPointCommunication() {
  preciceTrace1("~PointToPointCommunication()", _isConnected);

  closeConnection();
}

bool
PointToPointCommunication::isConnected() {
  return _isConnected;
}

void
PointToPointCommunication::acceptConnection(std::string const& nameAcceptor,
                                            std::string const& nameRequester) {
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);

  assertion(not isConnected());

  // Local (for process rank in the current participant) communication map that
  // defines a mapping from a process rank in the remote participant to an array
  // of local data indices, which define a subset of local (for process rank in
  // the current participant) data to be communicated between the current
  // process rank and the remote process rank.
  //
  // Example. Assume that the current process rank is 3. Assume that its
  // `communicationMap' is
  //
  //   1 -> {1, 3}
  //   4 -> {0, 2}
  //
  // then it means that the current process (with rank 3)
  // - has to communicate (send/receive) data with local indices 1 and 3 with
  //   the remote process with rank 1;
  // - has to communicate (send/receive) data with local indices 0 and 2 with
  //   the remote process with rank 4.
  std::map<int, std::vector<int>> communicationMap;

  {
    Event e("acceptorPreparation");

    if (utils::MasterSlave::_masterMode) {
      // Establish connection between participants' master processes.
      auto c = _communicationFactory->newCommunication();

      c->acceptConnection(nameAcceptor, nameRequester, 0, 1);
      // -----------------------------------------------------------------------
      // Exchange ranks of participants' master processes.
      int requesterMasterRank;

      c->send(utils::MasterSlave::_masterRank, 0);
      c->receive(requesterMasterRank, 0);
      // -----------------------------------------------------------------------
      // Exchange vertex distributions.
      auto& vertexDistribution = _mesh->getVertexDistribution();
      std::map<int, std::vector<int>> requesterVertexDistribution;

      m2n::send(c, vertexDistribution, 0);
      m2n::receive(c, requesterVertexDistribution, 0);
      // -----------------------------------------------------------------------
      // Iteratively construct different instances of communication map (from
      // the two vertex distributions), each of which corresponds to a single
      // acceptor process (in the current participant), and scatter them across
      // their corresponding acceptor processes (including the master acceptor
      // process).
      m2n::scatter(utils::MasterSlave::_communication,
                   communicationMap,
                   vertexDistribution,
                   requesterVertexDistribution);
    } else {
      assertion(utils::MasterSlave::_slaveMode);

      // Receive the corresponding (to the current acceptor process)
      // communication map that is being scattered by the master acceptor
      // process.
      m2n::receive(utils::MasterSlave::_communication, communicationMap, 0);
    }
  }

// NOTE:
// Change 0 to 1 to print `communicationMap'.
#if 0
  print(communicationMap);
#endif

  {
    Event e("acceptConnection");

    if (communicationMap.size() == 0) {
      _isConnected = true;

      return;
    }

    // Accept point-to-point connections between the current acceptor process
    // (in the current participant) with rank `utils::MasterSlave::_rank' and
    // (multiple) requester proccesses (in the requester participant).
    auto c = _communicationFactory->newCommunication();

    c->acceptConnectionAsServer(
        nameAcceptor + std::to_string(utils::MasterSlave::_rank),
        nameRequester,
        communicationMap.size());

    assertion(c->getRemoteCommunicatorSize() == communicationMap.size());

    _mappings.reserve(communicationMap.size());

    for (int localRequesterRank = 0;
         localRequesterRank < communicationMap.size();
         ++localRequesterRank) {
      int globalRequesterRank = -1;

      c->receive(globalRequesterRank, localRequesterRank);

      auto indices = std::move(communicationMap[globalRequesterRank]);

      // NOTE:
      // Everything is moved (efficiency)!
      _mappings.push_back(
          {localRequesterRank,
           globalRequesterRank,
           std::move(indices),
           // NOTE:
           // On the acceptor participant side, the communication object `c'
           // behaves as a server, i.e. it implicitly accepts multiple
           // connections to requester processes (in the requester participant)
           // according to point-to-point communicator size (in
           // `communicatorSizes') used during point-to-point connection
           // requests made by requester processes (in the requester
           // participant). As a result, only one communication object `c' is
           // needed to satisfy `communicationMap', and, therefore, for data
           // structure consistency of `_communications' with the requester
           // participant side, we simply duplicate references to the same
           // communication object `c'.
           c});
    }
  }

  _isConnected = true;
}

void
PointToPointCommunication::requestConnection(std::string const& nameAcceptor,
                                             std::string const& nameRequester) {
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);

  assertion(not isConnected());

  // Local (for process rank in the current participant) communication map that
  // defines a mapping from a process rank in the remote participant to an array
  // of local data indices, which define a subset of local (for process rank in
  // the current participant) data to be communicated between the current
  // process rank and the remote process rank.
  //
  // Example. Assume that the current process rank is 3. Assume that its
  // `communicationMap' is
  //
  //   1 -> {1, 3}
  //   4 -> {0, 2}
  //
  // then it means that the current process (with rank 3)
  // - has to communicate (send/receive) data with local indices 1 and 3 with
  //   the remote process with rank 1;
  // - has to communicate (send/receive) data with local indices 0 and 2 with
  //   the remote process with rank 4.
  std::map<int, std::vector<int>> communicationMap;

  {
    Event e("requesterPreparation");

    if (utils::MasterSlave::_masterMode) {
      // Establish connection between participants' master processes.
      auto c = _communicationFactory->newCommunication();

      c->requestConnection(nameAcceptor, nameRequester, 0, 1);
      // -----------------------------------------------------------------------
      // Exchange ranks of participants' master processes.
      int acceptorMasterRank;

      c->receive(acceptorMasterRank, 0);
      c->send(utils::MasterSlave::_masterRank, 0);
      // -----------------------------------------------------------------------
      // Exchange vertex distributions.
      auto& vertexDistribution = _mesh->getVertexDistribution();
      std::map<int, std::vector<int>> acceptorVertexDistribution;

      m2n::receive(c, acceptorVertexDistribution, 0);
      m2n::send(c, vertexDistribution, 0);
      // -----------------------------------------------------------------------
      // Iteratively construct different instances of communication map (from
      // the two vertex distributions), each of which corresponds to a single
      // requester process (in the current participant), and scatter them across
      // their corresponding requester processes (including the master requester
      // process).
      m2n::scatter(utils::MasterSlave::_communication,
                   communicationMap,
                   vertexDistribution,
                   acceptorVertexDistribution);
    } else {
      assertion(utils::MasterSlave::_slaveMode);

      // Receive the corresponding (to the current requester process)
      // communication map that is being scattered by the master requester
      // process.
      m2n::receive(utils::MasterSlave::_communication, communicationMap, 0);
    }
  }

// NOTE:
// Change 0 to 1 to print `communicationMap'.
#if 0
  print(communicationMap);
#endif

  {
    Event e("requestConnection");

    std::vector<com::PtrRequest> requests;

    requests.reserve(communicationMap.size());

    _mappings.reserve(communicationMap.size());



    // Request point-to-point connections between the current requester process
    // (in the current participant) and (multiple) acceptor proccesses (in the
    // acceptor participant) with ranks `acceptorRank' according to
    // communication map.
    for (auto& i : communicationMap) {
      auto globalAcceptorRank = i.first;
      auto indices = std::move(i.second);

      auto c = _communicationFactory->newCommunication();

      c->requestConnectionAsClient(
          nameAcceptor + std::to_string(globalAcceptorRank), nameRequester);

      assertion(c->getRemoteCommunicatorSize() == 1);

      auto request = c->aSend(&utils::MasterSlave::_rank, 0);

      requests.push_back(request);

      // NOTE:
      // Everything is moved (efficiency)!
      _mappings.push_back(
          {0,
           globalAcceptorRank,
           std::move(indices),
           // NOTE:
           // On the requester participant side, the communication objects
           // behave as clients, i.e. each of them requests only one connection
           // to acceptor process (in the acceptor participant).
           c});
    }

    for (auto request : requests) {
      request->wait();
    }
  }

  _isConnected = true;
}

void
PointToPointCommunication::closeConnection() {
  preciceTrace("closeConnection()");

  if (not isConnected())
    return;

  for (auto& mapping : _mappings) {
    mapping.communication->closeConnection();
  }

  _mappings.clear();

  _isConnected = false;
}

void
PointToPointCommunication::send(double* itemsToSend,
                                int size,
                                int valueDimension) {
  if (_mappings.size() == 0) {
    preciceCheck(
        size == 0, "send()", "Can't send anything from disconnected process!");

    return;
  }

  std::vector<com::PtrRequest> requests;

  requests.reserve(_mappings.size());

  std::vector<double> buffer;

  buffer.reserve(size * valueDimension);

  for (auto const& mapping : _mappings) {
    auto offset = buffer.size();

    for (auto index : mapping.indices) {
      for (int d = 0; d < valueDimension; ++d) {
        buffer.push_back(itemsToSend[index * valueDimension + d]);
      }
    }

    auto request =
        mapping.communication->aSend(buffer.data() + offset,
                                     mapping.indices.size() * valueDimension,
                                     mapping.localRemoteRank);

    requests.push_back(request);
  }

  for (auto request : requests) {
    request->wait();
  }
}

void
PointToPointCommunication::receive(double* itemsToReceive,
                                   int size,
                                   int valueDimension) {
  if (_mappings.size() == 0) {
    preciceCheck(size == 0,
                 "receive()",
                 "Can't receive anything to disconnected process!");

    return;
  }

  std::fill(itemsToReceive, itemsToReceive + size, 0);

  std::vector<double> buffer;

  buffer.reserve(size * valueDimension);

  for (auto const& mapping : _mappings) {
    auto offset = buffer.size();

    buffer.resize(buffer.size() + mapping.indices.size() * valueDimension);

    mapping.communication->receive(buffer.data() + offset,
                                   mapping.indices.size() * valueDimension,
                                   mapping.localRemoteRank);

    {
      int i = 0;

      for (auto index : mapping.indices) {
        for (int d = 0; d < valueDimension; ++d) {
          itemsToReceive[index * valueDimension + d] +=
              buffer[offset + i * valueDimension + d];
        }

        i++;
      }
    }
  }
}
}
} // namespace precice, m2n
