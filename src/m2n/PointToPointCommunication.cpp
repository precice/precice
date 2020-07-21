#include "PointToPointCommunication.hpp"
#include <algorithm>
#include <boost/container/flat_map.hpp>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <thread>
#include <vector>
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "com/CommunicationFactory.hpp"
#include "com/Request.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Event.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

using precice::utils::Event;

namespace precice {
bool extern syncMode;
namespace m2n {

void send(mesh::Mesh::VertexDistribution const &m,
          int                                   rankReceiver,
          com::PtrCommunication                 communication)
{
  communication->send(static_cast<int>(m.size()), rankReceiver);

  for (auto const &i : m) {
    auto const &rank    = i.first;
    auto const &indices = i.second;
    communication->send(rank, rankReceiver);
    communication->send(indices, rankReceiver);
  }
}

void receive(mesh::Mesh::VertexDistribution &m,
             int                             rankSender,
             com::PtrCommunication           communication)
{
  m.clear();
  int size = 0;
  communication->receive(size, rankSender);

  while (size--) {
    int rank = -1;
    communication->receive(rank, rankSender);
    communication->receive(m[rank], rankSender);
  }
}

void broadcastSend(mesh::Mesh::VertexDistribution const &m,
                   com::PtrCommunication                 communication = utils::MasterSlave::_communication)
{
  communication->broadcast(static_cast<int>(m.size()));

  for (auto const &i : m) {
    auto const &rank    = i.first;
    auto const &indices = i.second;
    communication->broadcast(rank);
    communication->broadcast(indices);
  }
}

void broadcastReceive(mesh::Mesh::VertexDistribution &m,
                      int                             rankBroadcaster,
                      com::PtrCommunication           communication = utils::MasterSlave::_communication)
{
  m.clear();
  int size = 0;
  communication->broadcast(size, rankBroadcaster);

  while (size--) {
    int rank = -1;
    communication->broadcast(rank, rankBroadcaster);
    communication->broadcast(m[rank], rankBroadcaster);
  }
}

void broadcast(mesh::Mesh::VertexDistribution &m)
{
  if (utils::MasterSlave::isMaster()) {
    // Broadcast (send) vertex distributions.
    m2n::broadcastSend(m);
  } else if (utils::MasterSlave::isSlave()) {
    // Broadcast (receive) vertex distributions.
    m2n::broadcastReceive(m, 0);
  }
}

void print(std::map<int, std::vector<int>> const &m)
{
  std::ostringstream oss;

  oss << "rank: " << utils::MasterSlave::getRank() << "\n";

  for (auto &i : m) {
    for (auto &j : i.second) {
      oss << i.first << ":" << j << '\n'; // prints rank:index
    }
  }

  if (utils::MasterSlave::isSlave()) {
    utils::MasterSlave::_communication->send(oss.str(), 0);
  } else {

    std::string s;

    for (int rank = 1; rank < utils::MasterSlave::getSize(); ++rank) {
      utils::MasterSlave::_communication->receive(s, rank);

      oss << s;
    }

    std::cout << oss.str();
  }
}

void printCommunicationPartnerCountStats(std::map<int, std::vector<int>> const &m)
{
  int size = m.size();

  if (utils::MasterSlave::isMaster()) {
    size_t count   = 0;
    size_t maximum = std::numeric_limits<size_t>::min();
    size_t minimum = std::numeric_limits<size_t>::max();
    size_t total   = size;

    if (size) {
      maximum = std::max(maximum, static_cast<size_t>(size));
      minimum = std::min(minimum, static_cast<size_t>(size));
      count++;
    }

    for (int rank = 1; rank < utils::MasterSlave::getSize(); ++rank) {
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

    auto average = static_cast<double>(total);
    if (count != 0) {
      average /= count;
    }

    std::cout << std::fixed << std::setprecision(3) //
              << "Number of Communication Partners per Interface Process:"
              << "\n"
              << "  Total:   " << total << "\n"
              << "  Maximum: " << maximum << "\n"
              << "  Minimum: " << minimum << "\n"
              << "  Average: " << average << "\n"
              << "Number of Interface Processes: " << count << "\n"
              << '\n';
  } else {
    PRECICE_ASSERT(utils::MasterSlave::isSlave());
    utils::MasterSlave::_communication->send(size, 0);
  }
}

void printLocalIndexCountStats(std::map<int, std::vector<int>> const &m)
{
  int size = 0;

  for (auto &i : m) {
    size += i.second.size();
  }

  if (utils::MasterSlave::isMaster()) {
    size_t count   = 0;
    size_t maximum = std::numeric_limits<size_t>::min();
    size_t minimum = std::numeric_limits<size_t>::max();
    size_t total   = size;

    if (size) {
      maximum = std::max(maximum, static_cast<size_t>(size));
      minimum = std::min(minimum, static_cast<size_t>(size));

      count++;
    }

    for (int rank = 1; rank < utils::MasterSlave::getSize(); ++rank) {
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

    auto average = static_cast<double>(total);
    if (count != 0) {
      average /= count;
    }

    std::cout << std::fixed << std::setprecision(3) //
              << "Number of LVDIs per Interface Process:"
              << "\n"
              << "  Total:   " << total << "\n"
              << "  Maximum: " << maximum << "\n"
              << "  Minimum: " << minimum << "\n"
              << "  Average: " << average << "\n"
              << "Number of Interface Processes: " << count << "\n"
              << '\n';
  } else {
    PRECICE_ASSERT(utils::MasterSlave::isSlave());

    utils::MasterSlave::_communication->send(size, 0);
  }
}

/** builds the communication map for a local distribution given the global distribution.
 *
 *
 * @param[in] thisVertexDistribution the local vertex distribution
 * @param[in] otherVertexDistribution the total vertex distribution
 * @param[in] thisRank the rank to build the map for
 *
 * @returns the resulting communication map for rank thisRank
 *
 * The approximate complexity of this function is:
 * \f$ \mathcal{O}(n \log(n) + m \log(n)) \f$
 *
 * * n is the total number of data indices for all ranks in `otherVertexDistribution'
 * * m is the number of local data indices for the current rank in `thisVertexDistribution`
 *
 */
std::map<int, std::vector<int>> buildCommunicationMap(
    // `thisVertexDistribution' is input vertex distribution from this participant.
    mesh::Mesh::VertexDistribution const &thisVertexDistribution,
    // `otherVertexDistribution' is input vertex distribution from other participant.
    mesh::Mesh::VertexDistribution const &otherVertexDistribution,
    int                                   thisRank = utils::MasterSlave::getRank())
{
  auto iterator = thisVertexDistribution.find(thisRank);
  if (iterator == thisVertexDistribution.end())
    return {};

  // Build lookup table from otherIndex -> rank for the otherVertexDistribution
  const auto lookupIndexRank = [&otherVertexDistribution] {
    boost::container::flat_multimap<int, int> lookupIndexRank;
    for (const auto &other : otherVertexDistribution) {
      for (const auto &otherIndex : other.second) {
        lookupIndexRank.emplace(otherIndex, other.first);
      }
    }
    return lookupIndexRank;
  }();

  auto const &indices = iterator->second;

  std::map<int, std::vector<int>> communicationMap;
  for (size_t index = 0lu; index < indices.size(); ++index) {
    auto range = lookupIndexRank.equal_range(indices[index]);
    for (auto iter = range.first; iter != range.second; ++iter) {
      communicationMap[iter->second].push_back(index);
    }
  }
  return communicationMap;
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
  PRECICE_TRACE(_isConnected);
  closeConnection();
}

bool PointToPointCommunication::isConnected() const
{
  return _isConnected;
}

void PointToPointCommunication::acceptConnection(std::string const &acceptorName,
                                                 std::string const &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected(), "Already connected.");

  mesh::Mesh::VertexDistribution &vertexDistribution = _mesh->getVertexDistribution();
  mesh::Mesh::VertexDistribution  requesterVertexDistribution;

  if (not utils::MasterSlave::isSlave()) {
    PRECICE_DEBUG("Exchange vertex distribution between both masters");
    Event e0("m2n.exchangeVertexDistribution");
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();

    c->acceptConnection(acceptorName, requesterName, "TMP-MASTERCOM-" + _mesh->getName(), utils::MasterSlave::getRank());

    // Exchange vertex distributions.
    m2n::send(vertexDistribution, 0, c);
    m2n::receive(requesterVertexDistribution, 0, c);
  }

  PRECICE_DEBUG("Broadcast vertex distributions");
  Event e1("m2n.broadcastVertexDistributions", precice::syncMode);
  m2n::broadcast(vertexDistribution);
  m2n::broadcast(requesterVertexDistribution);
  e1.stop();

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
  Event                           e2("m2n.buildCommunicationMap", precice::syncMode);
  std::map<int, std::vector<int>> communicationMap = m2n::buildCommunicationMap(
      vertexDistribution, requesterVertexDistribution);
  e2.stop();

// Print `communicationMap'.
#ifdef P2P_LCM_PRINT
  PRECICE_DEBUG("Print communication map");
  print(communicationMap);
#endif

// Print statistics of `communicationMap'.
#ifdef P2P_LCM_PRINT_STATS
  PRECICE_DEBUG("Print communication map statistics");
  printCommunicationPartnerCountStats(communicationMap);
  printLocalIndexCountStats(communicationMap);
#endif

  Event e4("m2n.createCommunications");
  e4.addData("Connections", communicationMap.size());
  if (communicationMap.empty()) {
    _isConnected = true;
    return;
  }

  PRECICE_DEBUG("Create and connect communication");
  _communication = _communicationFactory->newCommunication();

  // Accept point-to-point connections (as server) between the current acceptor
  // process (in the current participant) with rank `utils::MasterSlave::getRank()'
  // and (multiple) requester processes (in the requester participant).
  _communication->acceptConnectionAsServer(acceptorName,
                                           requesterName,
                                           _mesh->getName(),
                                           utils::MasterSlave::getRank(),
                                           communicationMap.size());

  PRECICE_DEBUG("Store communication map");
  for (auto const &comMap : communicationMap) {
    int  globalRequesterRank = comMap.first;
    auto indices             = std::move(communicationMap[globalRequesterRank]);

    _mappings.push_back({globalRequesterRank, std::move(indices), com::PtrRequest(), {}});
  }
  e4.stop();
  _isConnected = true;
}

void PointToPointCommunication::acceptPreConnection(std::string const &acceptorName,
                                                    std::string const &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected(), "Already connected.");

  const std::vector<int> &localConnectedRanks = _mesh->getConnectedRanks();

  if (localConnectedRanks.empty()) {
    _isConnected = true;
    return;
  }

  _communication = _communicationFactory->newCommunication();

  _communication->acceptConnectionAsServer(
      acceptorName,
      requesterName,
      _mesh->getName(),
      utils::MasterSlave::getRank(),
      localConnectedRanks.size());

  _connectionDataVector.reserve(localConnectedRanks.size());

  for (int connectedRank : localConnectedRanks) {
    _connectionDataVector.push_back({connectedRank, com::PtrRequest()});
  }

  _isConnected = true;
}

void PointToPointCommunication::requestConnection(std::string const &acceptorName,
                                                  std::string const &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected(), "Already connected.");

  mesh::Mesh::VertexDistribution &vertexDistribution = _mesh->getVertexDistribution();
  mesh::Mesh::VertexDistribution  acceptorVertexDistribution;

  if (not utils::MasterSlave::isSlave()) {
    PRECICE_DEBUG("Exchange vertex distribution between both masters");
    Event e0("m2n.exchangeVertexDistribution");
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();
    c->requestConnection(acceptorName, requesterName,
                         "TMP-MASTERCOM-" + _mesh->getName(),
                         0, 1);

    // Exchange vertex distributions.
    m2n::receive(acceptorVertexDistribution, 0, c);
    m2n::send(vertexDistribution, 0, c);
  }

  PRECICE_DEBUG("Broadcast vertex distributions");
  Event e1("m2n.broadcastVertexDistributions", precice::syncMode);
  m2n::broadcast(vertexDistribution);
  m2n::broadcast(acceptorVertexDistribution);
  e1.stop();

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
  Event                           e2("m2n.buildCommunicationMap", precice::syncMode);
  std::map<int, std::vector<int>> communicationMap = m2n::buildCommunicationMap(
      vertexDistribution, acceptorVertexDistribution);
  e2.stop();

// Print `communicationMap'.
#ifdef P2P_LCM_PRINT
  PRECICE_DEBUG("Print communication map");
  print(communicationMap);
#endif

// Print statistics of `communicationMap'.
#ifdef P2P_LCM_PRINT_STATS
  PRECICE_DEBUG("Print communication map statistics");
  printCommunicationPartnerCountStats(communicationMap);
  printLocalIndexCountStats(communicationMap);
#endif

  Event e4("m2n.createCommunications");
  e4.addData("Connections", communicationMap.size());
  if (communicationMap.empty()) {
    _isConnected = true;
    return;
  }

  std::vector<com::PtrRequest> requests;
  requests.reserve(communicationMap.size());

  std::set<int> acceptingRanks;
  for (auto &i : communicationMap)
    acceptingRanks.emplace(i.first);

  PRECICE_DEBUG("Create and connect communication");
  _communication = _communicationFactory->newCommunication();
  // Request point-to-point connections (as client) between the current
  // requester process (in the current participant) and (multiple) acceptor
  // processes (in the acceptor participant) to ranks `accceptingRanks'
  // according to `communicationMap`.
  _communication->requestConnectionAsClient(acceptorName, requesterName,
                                            _mesh->getName(),
                                            acceptingRanks, utils::MasterSlave::getRank());

  PRECICE_DEBUG("Store communication map");
  for (auto &i : communicationMap) {
    auto globalAcceptorRank = i.first;
    auto indices            = std::move(i.second);

    _mappings.push_back({globalAcceptorRank, std::move(indices), com::PtrRequest(), {}});
  }
  e4.stop();
  _isConnected = true;
}

void PointToPointCommunication::requestPreConnection(std::string const &acceptorName,
                                                     std::string const &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected(), "Already connected.");

  std::vector<int> localConnectedRanks = _mesh->getConnectedRanks();

  if (localConnectedRanks.empty()) {
    _isConnected = true;
    return;
  }

  std::vector<com::PtrRequest> requests;
  requests.reserve(localConnectedRanks.size());
  _connectionDataVector.reserve(localConnectedRanks.size());

  std::set<int> acceptingRanks(localConnectedRanks.begin(), localConnectedRanks.end());

  _communication = _communicationFactory->newCommunication();
  _communication->requestConnectionAsClient(acceptorName, requesterName,
                                            _mesh->getName(),
                                            acceptingRanks, utils::MasterSlave::getRank());

  for (auto &connectedRank : localConnectedRanks) {
    _connectionDataVector.push_back({connectedRank, com::PtrRequest()});
  }
  _isConnected = true;
}

void PointToPointCommunication::completeSlavesConnection()
{
  mesh::Mesh::CommunicationMap localCommunicationMap = _mesh->getCommunicationMap();

  for (auto &i : _connectionDataVector) {
    _mappings.push_back({i.remoteRank, std::move(localCommunicationMap[i.remoteRank]), i.request, {}});
  }
}

void PointToPointCommunication::closeConnection()
{
  PRECICE_TRACE();

  if (not isConnected())
    return;

  checkBufferedRequests(true);

  _communication.reset();
  _mappings.clear();
  _isConnected = false;
}

void PointToPointCommunication::send(double const *itemsToSend,
                                     size_t        size,
                                     int           valueDimension)
{

  if (_mappings.empty()) {
    return;
  }

  for (auto &mapping : _mappings) {
    // if (utils::MasterSlave::isMaster())
    //   std::cout<< "indices " << mapping.indices << std::endl;
    auto buffer = std::make_shared<std::vector<double>>();
    buffer->reserve(mapping.indices.size() * valueDimension);
    for (auto index : mapping.indices) {
      for (int d = 0; d < valueDimension; ++d) {
        buffer->push_back(itemsToSend[index * valueDimension + d]);
      }
    }
    auto request = _communication->aSend(*buffer, mapping.remoteRank);
    bufferedRequests.emplace_back(request, buffer);
  }
  checkBufferedRequests(false);
}

void PointToPointCommunication::receive(double *itemsToReceive,
                                        size_t  size,
                                        int     valueDimension)
{
  if (_mappings.empty()) {
    return;
  }

  std::fill(itemsToReceive, itemsToReceive + size, 0);

  for (auto &mapping : _mappings) {
    // if (not utils::MasterSlave::isMaster())
    //   std::cout<< "indices " << mapping.indices << std::endl;
    mapping.recvBuffer.resize(mapping.indices.size() * valueDimension);
    mapping.request = _communication->aReceive(mapping.recvBuffer, mapping.remoteRank);
  }

  for (auto &mapping : _mappings) {
    mapping.request->wait();

    int i = 0;
    for (auto index : mapping.indices) {
      for (int d = 0; d < valueDimension; ++d) {
        itemsToReceive[index * valueDimension + d] += mapping.recvBuffer[i * valueDimension + d];
      }
      i++;
    }
  }
}

void PointToPointCommunication::broadcastSend(const int &itemToSend)
{
  for (auto &connectionData : _connectionDataVector) {
    _communication->send(itemToSend, connectionData.remoteRank);
  }
}

void PointToPointCommunication::broadcastReceiveAll(std::vector<int> &itemToReceive)

{
  int data = 0;
  for (auto &connectionData : _connectionDataVector) {
    _communication->receive(data, connectionData.remoteRank);
    itemToReceive.push_back(data);
  }
}

void PointToPointCommunication::broadcastSendMesh()
{
  for (auto &connectionData : _connectionDataVector) {
    com::CommunicateMesh(_communication).sendMesh(*_mesh, connectionData.remoteRank);
  }
}

void PointToPointCommunication::broadcastReceiveAllMesh()
{
  for (auto &connectionData : _connectionDataVector) {
    com::CommunicateMesh(_communication).receiveMesh(*_mesh, connectionData.remoteRank);
  }
}

void PointToPointCommunication::scatterAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  for (auto &connectionData : _connectionDataVector) {
    _communication->send(localCommunicationMap[connectionData.remoteRank], connectionData.remoteRank);
  }
}

void PointToPointCommunication::gatherAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  for (auto &connectionData : _connectionDataVector) {
    _communication->receive(localCommunicationMap[connectionData.remoteRank], connectionData.remoteRank);
  }
}

void PointToPointCommunication::checkBufferedRequests(bool blocking)
{
  PRECICE_TRACE(bufferedRequests.size());
  do {
    for (auto it = bufferedRequests.begin(); it != bufferedRequests.end();) {
      if (it->first->test())
        it = bufferedRequests.erase(it);
      else
        ++it;
    }
    if (bufferedRequests.empty())
      return;
    if (blocking)
      std::this_thread::yield(); // give up our time slice, so MPI may work
  } while (blocking);
}

} // namespace m2n
} // namespace precice
