// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "PointToPointCommunication.hpp"

#include "mesh/Mesh.hpp"
#include "utils/EventTimings.hpp"
#include "utils/Globals.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Publisher.hpp"

#include <vector>

using precice::utils::Event;
using precice::utils::Publisher;

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
      if (not function(i.first, j))
        break;
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
send(std::vector<int> const& v,
     int rankReceiver,
     com::Communication::SharedPointer communication) {
  communication->send(static_cast<int>(v.size()), rankReceiver);
  communication->send(const_cast<int*>(&v[0]), v.size(), rankReceiver);
}

void
receive(std::vector<int>& v,
        int rankSender,
        com::Communication::SharedPointer communication) {
  v.clear();

  int size = 0;

  communication->receive(size, rankSender);

  v.resize(size);

  communication->receive(&v[0], size, rankSender);
}

void
send(std::map<int, std::vector<int>> const& m,
     int rankReceiver,
     com::Communication::SharedPointer communication) {
  communication->send(static_cast<int>(m.size()), rankReceiver);

  for (auto const& i : m) {
    auto const& rank = i.first;
    auto const& indices = i.second;

    communication->send(rank, rankReceiver);
    send(indices, rankReceiver, communication);
  }
}

void
receive(std::map<int, std::vector<int>>& m,
        int rankSender,
        com::Communication::SharedPointer communication) {
  m.clear();

  int size = 0;

  communication->receive(size, rankSender);

  while (size--) {
    int rank = -1;

    communication->receive(rank, rankSender);
    receive(m[rank], rankSender, communication);
  }
}

void
broadcast(std::vector<int> const& v,
          com::Communication::SharedPointer communication =
              utils::MasterSlave::_communication) {
  communication->broadcast(static_cast<int>(v.size()));
  communication->broadcast(const_cast<int*>(&v[0]), v.size());
}

void
broadcast(std::vector<int>& v,
          int rankBroadcaster,
          com::Communication::SharedPointer communication =
              utils::MasterSlave::_communication) {
  v.clear();

  int size = 0;

  communication->broadcast(size, rankBroadcaster);

  v.resize(size);

  communication->broadcast(&v[0], size, rankBroadcaster);
}

void
broadcastSend(std::map<int, std::vector<int>> const& m,
              com::Communication::SharedPointer communication =
                  utils::MasterSlave::_communication) {
  communication->broadcast(static_cast<int>(m.size()));

  for (auto const& i : m) {
    auto const& rank = i.first;
    auto const& indices = i.second;

    communication->broadcast(rank);
    broadcast(indices, communication);
  }
}

void
broadcastReceive(std::map<int, std::vector<int>>& m,
                 int rankBroadcaster,
                 com::Communication::SharedPointer communication =
                     utils::MasterSlave::_communication) {
  m.clear();

  int size = 0;

  communication->broadcast(size, rankBroadcaster);

  while (size--) {
    int rank = -1;

    communication->broadcast(rank, rankBroadcaster);
    broadcast(m[rank], rankBroadcaster, communication);
  }
}

void
broadcast(std::map<int, std::vector<int>>& m) {
  Event e(PointToPointCommunication::eventNamePrefix() +
              "PointToPointCommunication::broadcast",
          true);

  if (utils::MasterSlave::_masterMode) {
    // Broadcast (send) vertex distributions.
    m2n::broadcastSend(m);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    // Broadcast (receive) vertex distributions.
    m2n::broadcastReceive(m, 0);
  }
}

void
print(std::map<int, std::vector<int>> const& m) {
  std::ostringstream oss;

  oss << "rank: " << utils::MasterSlave::_rank << "\n";

  forMapOfRanges(m, [&oss](int rank, int index) {
    oss << rank << ":"
        << " " << index << "\n";

    return true;
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

// The approximate complexity of this function is O((number of local data
// indices for the current rank in `thisVertexDistribution') * (total number of
// data indices for all ranks in `otherVertexDistribution')).
std::map<int, std::vector<int>>
buildCommunicationMap(
    // `localIndexCount' is the number of unique local indices for the current
    // rank.
    size_t& localIndexCount,
    // `thisVertexDistribution' is input vertex distribution from this
    // participant.
    std::map<int, std::vector<int>> const& thisVertexDistribution,
    // `otherVertexDistribution' is input vertex distribution from other
    // participant.
    std::map<int, std::vector<int>> const& otherVertexDistribution,
    int thisRank = utils::MasterSlave::_rank) {
  Event e(PointToPointCommunication::eventNamePrefix() +
              "PointToPointCommunication::buildCommunicationMap",
          true);

  localIndexCount = 0;

  std::map<int, std::vector<int>> communicationMap;

  auto iterator = thisVertexDistribution.find(thisRank);

  if (iterator == thisVertexDistribution.end())
    return communicationMap;

  auto const& indices = iterator->second;

  localIndexCount = indices.size();

  int index = 0;

  forRange(indices, [&](int thisIndex) mutable {
    forMapOfRanges(
        otherVertexDistribution,
        [=, &communicationMap](int otherRank, int otherIndex) mutable {
          if (thisIndex == otherIndex) {
            communicationMap[otherRank].push_back(index);

            return false;
          }

          return true;
        });

    index++;
  });

  return communicationMap;
}

std::string PointToPointCommunication::_prefix;

PointToPointCommunication::ScopedSetEventNamePrefix::ScopedSetEventNamePrefix(
    std::string const& prefix)
    : _prefix(PointToPointCommunication::eventNamePrefix()) {
  PointToPointCommunication::setEventNamePrefix(prefix);
}

PointToPointCommunication::ScopedSetEventNamePrefix::
    ~ScopedSetEventNamePrefix() {
  PointToPointCommunication::setEventNamePrefix(_prefix);
}

void
PointToPointCommunication::setEventNamePrefix(std::string const& prefix) {
  _prefix = prefix;
}

std::string const&
PointToPointCommunication::eventNamePrefix() {
  return _prefix;
}

tarch::logging::Log PointToPointCommunication::_log(
    "precice::m2n::PointToPointCommunication");

PointToPointCommunication::PointToPointCommunication(
    com::CommunicationFactory::SharedPointer communicationFactory,
    mesh::PtrMesh mesh)
    : DistributedCommunication(mesh)
    , _communicationFactory(communicationFactory)
    , _localIndexCount(0)
    , _totalIndexCount(0)
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

  preciceCheck(not isConnected(), "acceptConnection()", "Already connected!");

  Event e(_prefix + "PointToPointCommunication::acceptConnection", true);

  std::map<int, std::vector<int>>& vertexDistribution =
      _mesh->getVertexDistribution();
  std::map<int, std::vector<int>> requesterVertexDistribution;

  if (utils::MasterSlave::_masterMode) {
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();

    {
      Event e(
          _prefix + "PointToPointCommunication::acceptConnection/synchronize",
          true);

      c->acceptConnection(nameAcceptor, nameRequester, 0, 1);
    }

    Event e(_prefix + "PointToPointCommunication::acceptConnection/exchange",
            true);

    int requesterMasterRank;

    // Exchange ranks of participants' master processes.
    c->send(utils::MasterSlave::_masterRank, 0);
    c->receive(requesterMasterRank, 0);

    // Exchange vertex distributions.
    m2n::send(vertexDistribution, 0, c);
    m2n::receive(requesterVertexDistribution, 0, c);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    Event(_prefix + "PointToPointCommunication::acceptConnection/synchronize",
          true);
    Event(_prefix + "PointToPointCommunication::acceptConnection/exchange",
          true);
  }

  m2n::broadcast(vertexDistribution);
  m2n::broadcast(requesterVertexDistribution);

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
  std::map<int, std::vector<int>> communicationMap = m2n::buildCommunicationMap(
      _localIndexCount, vertexDistribution, requesterVertexDistribution);

// NOTE:
// Change 0 to 1 to print `communicationMap'.
#if 0
  print(communicationMap);
#endif

#ifdef SuperMUC_WORK
  try {
    auto addressDirectory = _communicationFactory->addressDirectory();

    if (utils::MasterSlave::_masterMode) {
      Event e(_prefix +
              "PointToPointCommunication::acceptConnection/createDirectories");

      for (int rank = 0; rank < utils::MasterSlave::_size; ++rank) {
        Publisher::createDirectory(addressDirectory + "/" + "." + nameAcceptor +
                                   "-" + std::to_string(rank) + ".address");
      }
    }

    utils::Parallel::synchronizeProcesses();
  } catch (...) {
  }
#endif

  Event e2(_prefix + "PointToPointCommunication::acceptConnection/accept",
           true);

  if (communicationMap.size() == 0) {
    _isConnected = true;

    return;
  }

  // Accept point-to-point connections (as server) between the current acceptor
  // process (in the current participant) with rank `utils::MasterSlave::_rank'
  // and (multiple) requester proccesses (in the requester participant).
  auto c = _communicationFactory->newCommunication();

#ifdef SuperMUC_WORK
  Publisher::ScopedPushDirectory spd("." + nameAcceptor + "-" +
                                     std::to_string(utils::MasterSlave::_rank) +
                                     ".address");
#endif

  c->acceptConnectionAsServer(
      nameAcceptor + "-" + std::to_string(utils::MasterSlave::_rank),
      nameRequester,
      communicationMap.size());

  assertion(c->getRemoteCommunicatorSize() == communicationMap.size());

  _mappings.reserve(communicationMap.size());

  for (int localRequesterRank = 0; localRequesterRank < communicationMap.size();
       ++localRequesterRank) {
    int globalRequesterRank = -1;

    c->receive(globalRequesterRank, localRequesterRank);

    auto indices = std::move(communicationMap[globalRequesterRank]);

    _totalIndexCount += indices.size();

    // NOTE:
    // Everything is moved (efficiency)!
    _mappings.push_back(
        {localRequesterRank,
         globalRequesterRank,
         std::move(indices),
         // NOTE:
         // On the acceptor participant side, the communication object `c'
         // behaves as a server, i.e. it implicitly accepts multiple connections
         // to requester processes (in the requester participant). As a result,
         // only one communication object `c' is needed to satisfy
         // `communicationMap', and, therefore, for data structure consistency
         // of `_mappings' with the requester participant side, we simply
         // duplicate references to the same communication object `c'.
         c});
  }

  _requests.reserve(_mappings.size());

  _buffer.reserve(_totalIndexCount * _mesh->getDimensions());

  _isConnected = true;
}

void
PointToPointCommunication::requestConnection(std::string const& nameAcceptor,
                                             std::string const& nameRequester) {
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);

  preciceCheck(not isConnected(), "requestConnection()", "Already connected!");

  Event e(_prefix + "PointToPointCommunication::requestConnection", true);

  std::map<int, std::vector<int>>& vertexDistribution =
      _mesh->getVertexDistribution();
  std::map<int, std::vector<int>> acceptorVertexDistribution;

  if (utils::MasterSlave::_masterMode) {
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();

    {
      Event e(
          _prefix + "PointToPointCommunication::requestConnection/synchronize",
          true);

      Publisher::ScopedSetEventNamePrefix ssenp(
          _prefix +
          "PointToPointCommunication::requestConnection"
          "/"
          "synchronize"
          "/");

      c->requestConnection(nameAcceptor, nameRequester, 0, 1);
    }

    Event e(_prefix + "PointToPointCommunication::requestConnection/exchange",
            true);

    int acceptorMasterRank;

    // Exchange ranks of participants' master processes.
    c->receive(acceptorMasterRank, 0);
    c->send(utils::MasterSlave::_masterRank, 0);

    // Exchange vertex distributions.
    m2n::receive(acceptorVertexDistribution, 0, c);
    m2n::send(vertexDistribution, 0, c);
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    Event(_prefix + "PointToPointCommunication::requestConnection/synchronize",
          true);
    Event(_prefix + "PointToPointCommunication::requestConnection/exchange",
          true);
  }

  m2n::broadcast(vertexDistribution);
  m2n::broadcast(acceptorVertexDistribution);

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
  std::map<int, std::vector<int>> communicationMap = m2n::buildCommunicationMap(
      _localIndexCount, vertexDistribution, acceptorVertexDistribution);

// NOTE:
// Change 0 to 1 to print `communicationMap'.
#if 0
  print(communicationMap);
#endif

#ifdef SuperMUC_WORK
  try {
    auto addressDirectory = _communicationFactory->addressDirectory();

    utils::Parallel::synchronizeProcesses();
  } catch (...) {
  }
#endif

  Event e2(_prefix + "PointToPointCommunication::requestConnection/request",
           true);

  Publisher::ScopedSetEventNamePrefix ssenp(
      _prefix +
      "PointToPointCommunication::requestConnection"
      "/"
      "request"
      "/");

  std::vector<com::Request::SharedPointer> requests;

  requests.reserve(communicationMap.size());

  _mappings.reserve(communicationMap.size());

  // Request point-to-point connections (as client) between the current
  // requester process (in the current participant) and (multiple) acceptor
  // proccesses (in the acceptor participant) with ranks `globalAcceptorRank'
  // according to communication map.
  for (auto& i : communicationMap) {
    auto globalAcceptorRank = i.first;
    auto indices = std::move(i.second);

    _totalIndexCount += indices.size();

    auto c = _communicationFactory->newCommunication();

#ifdef SuperMUC_WORK
    Publisher::ScopedPushDirectory spd("." + nameAcceptor + "-" +
                                       std::to_string(globalAcceptorRank) +
                                       ".address");
#endif

    c->requestConnectionAsClient(
        nameAcceptor + "-" + std::to_string(globalAcceptorRank), nameRequester);

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
         // On the requester participant side, the communication objects behave
         // as clients, i.e. each of them requests only one connection to
         // acceptor process (in the acceptor participant).
         c});
  }

  com::Request::wait(requests);

  _requests.reserve(_mappings.size());

  _buffer.reserve(_totalIndexCount * _mesh->getDimensions());

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

  _requests.clear();

  _buffer.clear();

  _localIndexCount = 0;

  _totalIndexCount = 0;

  _isConnected = false;
}

void
PointToPointCommunication::send(double* itemsToSend,
                                int size,
                                int valueDimension) {
  Event e(_prefix + "PointToPointCommunication::send", true);

  preciceInfo("send(double)",
              "Size"
                  << ":"
                  << " " << size << ";"
                  << " "
                  << "Dimension"
                  << ":"
                  << " " << valueDimension << ".");

  if (_mappings.size() == 0) {
    preciceCheck(size == 0 && _localIndexCount == 0,
                 "send()",
                 "Can't send anything from disconnected process!");

    return;
  }

  preciceCheck(size == _localIndexCount * valueDimension,
               "send()",
               "Inconsistency between expected"
                   << " "
                   << "(" << _localIndexCount * valueDimension << ")"
                   << " "
                   << "and provided"
                   << " "
                   << "(" << size << ")"
                   << " "
                   << "data sizes!");

  for (auto const& mapping : _mappings) {
    auto offset = _buffer.size();

    for (auto index : mapping.indices) {
      for (int d = 0; d < valueDimension; ++d) {
        _buffer.push_back(itemsToSend[index * valueDimension + d]);
      }
    }

    auto request =
        mapping.communication->aSend(_buffer.data() + offset,
                                     mapping.indices.size() * valueDimension,
                                     mapping.localRemoteRank);

    _requests.push_back(request);
  }

  com::Request::wait(_requests);

  _requests.clear();

  _buffer.clear();
}

void
PointToPointCommunication::receive(double* itemsToReceive,
                                   int size,
                                   int valueDimension) {
  Event e(_prefix + "PointToPointCommunication::receive", true);

  preciceInfo("receive(double)",
              "Size"
                  << ":"
                  << " " << size << ";"
                  << " "
                  << "Dimension"
                  << ":"
                  << " " << valueDimension << ".");

  if (_mappings.size() == 0) {
    preciceCheck(size == 0 && _localIndexCount == 0,
                 "receive()",
                 "Can't receive anything to disconnected process!");

    return;
  }

  preciceCheck(size == _localIndexCount * valueDimension,
               "receive()",
               "Inconsistency between expected"
                   << " "
                   << "(" << _localIndexCount * valueDimension << ")"
                   << " "
                   << "and provided"
                   << " "
                   << "(" << size << ")"
                   << " "
                   << "data sizes!");

  std::fill(itemsToReceive, itemsToReceive + size, 0);

  for (auto const& mapping : _mappings) {
    auto offset = _buffer.size();

    _buffer.resize(_buffer.size() + mapping.indices.size() * valueDimension);

    mapping.communication->receive(_buffer.data() + offset,
                                   mapping.indices.size() * valueDimension,
                                   mapping.localRemoteRank);

    {
      int i = 0;

      for (auto index : mapping.indices) {
        for (int d = 0; d < valueDimension; ++d) {
          itemsToReceive[index * valueDimension + d] +=
              _buffer[offset + i * valueDimension + d];
        }

        i++;
      }
    }
  }

  _requests.clear();

  _buffer.clear();
}
}
} // namespace precice, m2n
