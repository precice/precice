#include "PointToPointCommunication.hpp"
#include <vector>
#include "com/Communication.hpp"
#include "com/CommunicationFactory.hpp"
#include "mesh/Mesh.hpp"
#include "utils/EventTimings.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Publisher.hpp"

using precice::utils::Event;
using precice::utils::Publisher;

namespace precice
{
namespace m2n
{

void send(std::vector<int> const &v,
          int                     rankReceiver,
          com::PtrCommunication   communication)
{
  communication->send(static_cast<int>(v.size()), rankReceiver);
  communication->send(const_cast<int *>(&v[0]), v.size(), rankReceiver);
}

void receive(std::vector<int> &    v,
             int                   rankSender,
             com::PtrCommunication communication)
{
  v.clear();
  int size = 0;
  communication->receive(size, rankSender);
  v.resize(size);
  communication->receive(v.data(), size, rankSender);
}

void send(std::map<int, std::vector<int>> const &m,
          int                                    rankReceiver,
          com::PtrCommunication                  communication)
{
  communication->send(static_cast<int>(m.size()), rankReceiver);

  for (auto const &i : m) {
    auto const &rank    = i.first;
    auto const &indices = i.second;
    communication->send(rank, rankReceiver);
    send(indices, rankReceiver, communication);
  }
}

void receive(std::map<int, std::vector<int>> &m,
             int                              rankSender,
             com::PtrCommunication            communication)
{
  m.clear();
  int size = 0;
  communication->receive(size, rankSender);

  while (size--) {
    int rank = -1;
    communication->receive(rank, rankSender);
    receive(m[rank], rankSender, communication);
  }
}

void broadcast(std::vector<int> const &v,
               com::PtrCommunication   communication = utils::MasterSlave::_communication)
{
  communication->broadcast(static_cast<int>(v.size()));
  communication->broadcast(const_cast<int *>(&v[0]), v.size());
}

void broadcast(std::vector<int> &    v,
               int                   rankBroadcaster,
               com::PtrCommunication communication = utils::MasterSlave::_communication)
{
  v.clear();
  int size = 0;
  communication->broadcast(size, rankBroadcaster);
  v.resize(size);
  communication->broadcast(&v[0], size, rankBroadcaster);
}

void broadcastSend(std::map<int, std::vector<int>> const &m,
                   com::PtrCommunication                  communication = utils::MasterSlave::_communication)
{
  communication->broadcast(static_cast<int>(m.size()));

  for (auto const &i : m) {
    auto const &rank    = i.first;
    auto const &indices = i.second;
    communication->broadcast(rank);
    broadcast(indices, communication);
  }
}

void broadcastReceive(std::map<int, std::vector<int>> &m,
                      int                              rankBroadcaster,
                      com::PtrCommunication            communication = utils::MasterSlave::_communication)
{
  m.clear();
  int size = 0;
  communication->broadcast(size, rankBroadcaster);

  while (size--) {
    int rank = -1;
    communication->broadcast(rank, rankBroadcaster);
    broadcast(m[rank], rankBroadcaster, communication);
  }
}

void broadcast(std::map<int, std::vector<int>> &m)
{
  if (utils::MasterSlave::_masterMode) {
    // Broadcast (send) vertex distributions.
    m2n::broadcastSend(m);
  } else {
    assertion(utils::MasterSlave::_slaveMode);
    // Broadcast (receive) vertex distributions.
    m2n::broadcastReceive(m, 0);
  }
}

void print(std::map<int, std::vector<int>> const &m)
{
  std::ostringstream oss;

  oss << "rank: " << utils::MasterSlave::_rank << "\n";

  for (auto &i : m) {
    for (auto &j : i.second) {
      oss << i.first << ":" << j << std::endl; // prints rank:index
    }
  }

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

void printCommunicationPartnerCountStats(std::map<int, std::vector<int>> const &m)
{
  int size = m.size();

  if (utils::MasterSlave::_masterMode) {
    size_t count   = 0;
    size_t maximum = std::numeric_limits<size_t>::min();
    size_t minimum = std::numeric_limits<size_t>::max();
    size_t total   = size;

    if (size) {
      maximum = std::max(maximum, static_cast<size_t>(size));
      minimum = std::min(minimum, static_cast<size_t>(size));
      count++;
    }

    for (int rank = 1; rank < utils::MasterSlave::_size; ++rank) {
      utils::MasterSlave::_communication->receive(size, rank);

      total += size;

      if (size) {
        maximum = std::max(maximum, static_cast<size_t>(size));
        minimum = std::min(minimum, static_cast<size_t>(size));
        count++;
      }
    }

    if (minimum > maximum)
      minimum = maximum;

    auto average = static_cast<double>(total) / count;

    std::cout << std::fixed << std::setprecision(3) //
              << "Number of Communication Partners per Interface Process:"
              << "\n"
              << "  Total:   " << total << "\n"
              << "  Maximum: " << maximum << "\n"
              << "  Minimum: " << minimum << "\n"
              << "  Average: " << average << "\n"
              << "Number of Interface Processes: " << count << "\n"
              << std::endl;
  } else {
    assertion(utils::MasterSlave::_slaveMode);
    utils::MasterSlave::_communication->send(size, 0);
  }
}

void printLocalIndexCountStats(std::map<int, std::vector<int>> const &m)
{
  int size = 0;

  for (auto &i : m) {
    size += i.second.size();
  }

  if (utils::MasterSlave::_masterMode) {
    size_t count   = 0;
    size_t maximum = std::numeric_limits<size_t>::min();
    size_t minimum = std::numeric_limits<size_t>::max();
    size_t total   = size;

    if (size) {
      maximum = std::max(maximum, static_cast<size_t>(size));
      minimum = std::min(minimum, static_cast<size_t>(size));

      count++;
    }

    for (int rank = 1; rank < utils::MasterSlave::_size; ++rank) {
      utils::MasterSlave::_communication->receive(size, rank);

      total += size;

      if (size) {
        maximum = std::max(maximum, static_cast<size_t>(size));
        minimum = std::min(minimum, static_cast<size_t>(size));

        count++;
      }
    }

    if (minimum > maximum)
      minimum = maximum;

    auto average = static_cast<double>(total) / count;

    std::cout << std::fixed << std::setprecision(3) //
              << "Number of LVDIs per Interface Process:"
              << "\n"
              << "  Total:   " << total << "\n"
              << "  Maximum: " << maximum << "\n"
              << "  Minimum: " << minimum << "\n"
              << "  Average: " << average << "\n"
              << "Number of Interface Processes: " << count << "\n"
              << std::endl;
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    utils::MasterSlave::_communication->send(size, 0);
  }
}

// The approximate complexity of this function is O((number of local data
// indices for the current rank in `thisVertexDistribution') * (total number of
// data indices for all ranks in `otherVertexDistribution')).
std::map<int, std::vector<int>> buildCommunicationMap(
    // `localIndexCount' is the number of unique local indices for the current rank.
    size_t &localIndexCount,
    // `thisVertexDistribution' is input vertex distribution from this participant.
    std::map<int, std::vector<int>> const &thisVertexDistribution,
    // `otherVertexDistribution' is input vertex distribution from other participant.
    std::map<int, std::vector<int>> const &otherVertexDistribution,
    int                                    thisRank = utils::MasterSlave::_rank)
{

  localIndexCount = 0;

  std::map<int, std::vector<int>> communicationMap;

  auto iterator = thisVertexDistribution.find(thisRank);

  if (iterator == thisVertexDistribution.end())
    return communicationMap;

  auto const &indices = iterator->second;

  int index = 0;

  for (int thisIndex : indices) {
    for (const auto &other : otherVertexDistribution) {
      for (const auto &otherIndex : other.second) {
        if (thisIndex == otherIndex) {
          communicationMap[other.first].push_back(index);
          break;
        }
      }
    }
    ++index;
  }

  // CAUTION:
  // This prevents point-to-point communication from considering those process
  // ranks, which don't have matching indices in the remote participant
  // (i.e. are not supposed to have any communication partners at all). For
  // instance, this happens a lot in case of unstructured grids on the
  // interface, or because user did a mistake and did not define the interface
  // boundary properly.
  if (communicationMap.size() > 0)
    localIndexCount = indices.size();

  return communicationMap;
}

std::string PointToPointCommunication::_prefix;

PointToPointCommunication::ScopedSetEventNamePrefix::ScopedSetEventNamePrefix(
    std::string const &prefix)
    : _prefix(PointToPointCommunication::eventNamePrefix())
{
  PointToPointCommunication::setEventNamePrefix(prefix);
}

PointToPointCommunication::ScopedSetEventNamePrefix::~ScopedSetEventNamePrefix()
{
  PointToPointCommunication::setEventNamePrefix(_prefix);
}

void PointToPointCommunication::setEventNamePrefix(std::string const &prefix)
{
  _prefix = prefix;
}

std::string const &PointToPointCommunication::eventNamePrefix()
{
  return _prefix;
}


PointToPointCommunication::PointToPointCommunication(
    com::PtrCommunicationFactory communicationFactory,
    mesh::PtrMesh                mesh)
    : DistributedCommunication(mesh),
      _communicationFactory(communicationFactory)
{
}

PointToPointCommunication::~PointToPointCommunication()
{
  TRACE(_isConnected);
  closeConnection();
}

bool PointToPointCommunication::isConnected()
{
  return _isConnected;
}

void PointToPointCommunication::acceptConnection(std::string const &nameAcceptor,
                                                 std::string const &nameRequester)
{
  TRACE(nameAcceptor, nameRequester);
  CHECK(not isConnected(), "Already connected!");
  CHECK(utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode,
        "You can only use a point-to-point communication between two participants which both use a master. "
            << "Please use distribution-type gather-scatter instead.");

  std::map<int, std::vector<int>> &vertexDistribution = _mesh->getVertexDistribution();
  std::map<int, std::vector<int>>  requesterVertexDistribution;

  if (utils::MasterSlave::_masterMode) {
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();

    c->acceptConnection(nameAcceptor, nameRequester, 1);

    int requesterMasterRank;

    // Exchange ranks of participants' master processes.
    c->send(utils::MasterSlave::_masterRank, 0);
    c->receive(requesterMasterRank, 0);

    // Exchange vertex distributions.
    m2n::send(vertexDistribution, 0, c);
    m2n::receive(requesterVertexDistribution, 0, c);
  } else {
    assertion(utils::MasterSlave::_slaveMode);
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

// Print `communicationMap'.
#ifdef P2P_LCM_PRINT
  e.stop(true);
  print(communicationMap);
  e.start(true);
#endif

// Print statistics of `communicationMap'.
#ifdef P2P_LCM_PRINT_STATS
  e.stop(true);
  printCommunicationPartnerCountStats(communicationMap);
  printLocalIndexCountStats(communicationMap);
  e.start(true);
#endif

#ifdef SuperMUC_WORK
  try {
    auto addressDirectory = _communicationFactory->addressDirectory();

    if (utils::MasterSlave::_masterMode) {
      Event e(_prefix + "PointToPointCommunication::acceptConnection/createDirectories");

      for (int rank = 0; rank < utils::MasterSlave::_size; ++rank) {
        Publisher::createDirectory(addressDirectory + "/" + "." + nameAcceptor + "-" + _mesh->getName() +
                                   "-" + std::to_string(rank) + ".address");
      }
    }
    utils::Parallel::synchronizeProcesses();
  } catch (...) {
  }
#endif

  if (communicationMap.empty()) {
    assertion(_localIndexCount == 0);
    _isConnected = true;
    return;
  }

  // Accept point-to-point connections (as server) between the current acceptor
  // process (in the current participant) with rank `utils::MasterSlave::_rank'
  // and (multiple) requester processes (in the requester participant).
  auto c = _communicationFactory->newCommunication();

#ifdef SuperMUC_WORK
  Publisher::ScopedPushDirectory spd("." + nameAcceptor + "-" + _mesh->getName() + "-" +
                                     std::to_string(utils::MasterSlave::_rank) + ".address");
#endif

  c->acceptConnectionAsServer(
      nameAcceptor + "-" + std::to_string(utils::MasterSlave::_rank),
      nameRequester,
      communicationMap.size());

  // assertion(c->getRemoteCommunicatorSize() == communicationMap.size());

  _mappings.reserve(communicationMap.size());

  for (size_t localRequesterRank = 0; localRequesterRank < communicationMap.size(); ++localRequesterRank) {
    int globalRequesterRank = -1;

    c->receive(globalRequesterRank, localRequesterRank);

    auto indices = std::move(communicationMap[globalRequesterRank]);
    _totalIndexCount += indices.size();

    // NOTE:
    // Everything is moved (efficiency)!
    _mappings.push_back(
        {static_cast<int>(localRequesterRank),
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

  _buffer.reserve(_totalIndexCount * _mesh->getDimensions());
  _isConnected = true;
}

void PointToPointCommunication::requestConnection(std::string const &nameAcceptor,
                                                  std::string const &nameRequester)
{
  TRACE(nameAcceptor, nameRequester);
  CHECK(not isConnected(), "Already connected!");
  CHECK(utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode,
        "You can only use a point-to-point communication between two participants which both use a master. "
            << "Please use distribution-type gather-scatter instead.");

  std::map<int, std::vector<int>> &vertexDistribution = _mesh->getVertexDistribution();
  std::map<int, std::vector<int>>  acceptorVertexDistribution;

  if (utils::MasterSlave::_masterMode) {
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();
    {
      Publisher::ScopedSetEventNamePrefix ssenp(
          _prefix + "PointToPointCommunication::requestConnection/synchronize/");

      c->requestConnection(nameAcceptor, nameRequester, 0, 1);
    }

    int acceptorMasterRank;

    // Exchange ranks of participants' master processes.
    c->receive(acceptorMasterRank, 0);
    c->send(utils::MasterSlave::_masterRank, 0);

    // Exchange vertex distributions.
    m2n::receive(acceptorVertexDistribution, 0, c);
    m2n::send(vertexDistribution, 0, c);
  } else {
    assertion(utils::MasterSlave::_slaveMode);
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

// Print `communicationMap'.
#ifdef P2P_LCM_PRINT
  e.stop(true);
  print(communicationMap);
  e.start(true);
#endif

// Print statistics of `communicationMap'.
#ifdef P2P_LCM_PRINT_STATS
  e.stop(true);
  printCommunicationPartnerCountStats(communicationMap);
  printLocalIndexCountStats(communicationMap);
  e.start(true);
#endif

#ifdef SuperMUC_WORK
  try {
    auto addressDirectory = _communicationFactory->addressDirectory();

    utils::Parallel::synchronizeProcesses();
  } catch (...) {
  }
#endif

  if (communicationMap.empty()) {
    assertion(_localIndexCount == 0);
    _isConnected = true;
    return;
  }

  Publisher::ScopedSetEventNamePrefix ssenp(
      _prefix + "PointToPointCommunication::requestConnection/request/");

  std::vector<com::PtrRequest> requests;
  requests.reserve(communicationMap.size());
  _mappings.reserve(communicationMap.size());

  // Request point-to-point connections (as client) between the current
  // requester process (in the current participant) and (multiple) acceptor
  // processes (in the acceptor participant) with ranks `globalAcceptorRank'
  // according to communication map.
  for (auto &i : communicationMap) {
    auto globalAcceptorRank = i.first;
    auto indices            = std::move(i.second);

    _totalIndexCount += indices.size();

    auto c = _communicationFactory->newCommunication();

#ifdef SuperMUC_WORK
    Publisher::ScopedPushDirectory spd("." + nameAcceptor + "-" + _mesh->getName() + "-" +
                                       std::to_string(globalAcceptorRank) + ".address");
#endif

    c->requestConnectionAsClient(nameAcceptor + "-" + std::to_string(globalAcceptorRank), nameRequester);
    // assertion(c->getRemoteCommunicatorSize() == 1);

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
  _buffer.reserve(_totalIndexCount * _mesh->getDimensions());
  _isConnected = true;
}

void PointToPointCommunication::closeConnection()
{
  TRACE();

  if (not isConnected())
    return;

  for (auto &mapping : _mappings) {
    mapping.communication->closeConnection();
  }

  _mappings.clear();
  _buffer.clear();
  _localIndexCount = 0;
  _totalIndexCount = 0;
  _isConnected     = false;
}

void PointToPointCommunication::send(double *itemsToSend,
                                     size_t  size,
                                     int     valueDimension)
{

  if (_mappings.size() == 0) {
    assertion(_localIndexCount == 0);
    return;
  }

  assertion(size == _localIndexCount * valueDimension, size, _localIndexCount * valueDimension);

  for (auto &mapping : _mappings) {
    mapping.offset = _buffer.size();

    for (auto index : mapping.indices) {
      for (int d = 0; d < valueDimension; ++d) {
        _buffer.push_back(itemsToSend[index * valueDimension + d]);
      }
    }

    mapping.request = mapping.communication->aSend(_buffer.data() + mapping.offset,
                                                   mapping.indices.size() * valueDimension,
                                                   mapping.localRemoteRank);
  }

  for (auto &mapping : _mappings) {
    mapping.request->wait();
  }
  _buffer.clear();
}

void PointToPointCommunication::receive(double *itemsToReceive,
                                        size_t  size,
                                        int     valueDimension)
{
  if (_mappings.size() == 0) {
    assertion(_localIndexCount == 0);
    return;
  }
  assertion(size == _localIndexCount * valueDimension, size, _localIndexCount * valueDimension);

  std::fill(itemsToReceive, itemsToReceive + size, 0);

  for (auto &mapping : _mappings) {
    mapping.offset = _buffer.size();
    _buffer.resize(_buffer.size() + mapping.indices.size() * valueDimension);

    mapping.request =
        mapping.communication->aReceive(_buffer.data() + mapping.offset,
                                        mapping.indices.size() * valueDimension,
                                        mapping.localRemoteRank);
  }

  for (auto &mapping : _mappings) {
    mapping.request->wait();

    int i = 0;

    for (auto index : mapping.indices) {
      for (int d = 0; d < valueDimension; ++d) {
        itemsToReceive[index * valueDimension + d] += _buffer[mapping.offset + i * valueDimension + d];
      }

      i++;
    }
  }

  _buffer.clear();
}
} // namespace m2n
} // namespace precice
