#include "PointToPointCommunication.hpp"
#include <iomanip>
#include <vector>
#include <thread>
#include "com/Communication.hpp"
#include "com/CommunicationFactory.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Event.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Publisher.hpp"

using precice::utils::Event;
using precice::utils::Publisher;

namespace precice
{
bool extern syncMode;
namespace m2n
{

void send(mesh::Mesh::VertexDistribution const &m,
          int                                    rankReceiver,
          com::PtrCommunication                  communication)
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
             int                              rankSender,
             com::PtrCommunication            communication)
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
                   com::PtrCommunication                  communication = utils::MasterSlave::_communication)
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
                      int                              rankBroadcaster,
                      com::PtrCommunication            communication = utils::MasterSlave::_communication)
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
      oss << i.first << ":" << j << '\n'; // prints rank:index
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
              << '\n';
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
              << '\n';
  } else {
    assertion(utils::MasterSlave::_slaveMode);

    utils::MasterSlave::_communication->send(size, 0);
  }
}

// The approximate complexity of this function is O((number of local data
// indices for the current rank in `thisVertexDistribution') * (total number of
// data indices for all ranks in `otherVertexDistribution')).
std::map<int, std::vector<int>> buildCommunicationMap(
    // `thisVertexDistribution' is input vertex distribution from this participant.
    mesh::Mesh::VertexDistribution const &thisVertexDistribution,
    // `otherVertexDistribution' is input vertex distribution from other participant.
    mesh::Mesh::VertexDistribution const &otherVertexDistribution,
    int                                    thisRank = utils::MasterSlave::_rank)
{
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
  TRACE(_isConnected);
  closeConnection();
}

bool PointToPointCommunication::isConnected()
{
  return _isConnected;
}

void PointToPointCommunication::acceptConnection(std::string const &acceptorName,
                                                 std::string const &requesterName)
{
  TRACE(acceptorName, requesterName);
  CHECK(not isConnected(), "Already connected!");
  CHECK(utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode,
        "You can only use a point-to-point communication between two participants which both use a master. "
            << "Please use distribution-type gather-scatter instead.");

  mesh::Mesh::VertexDistribution &vertexDistribution = _mesh->getVertexDistribution();
  mesh::Mesh::VertexDistribution  requesterVertexDistribution;

  if (utils::MasterSlave::_masterMode) {
    Event e0("m2n.exchangeVertexDistribution");
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();

    c->acceptConnection(acceptorName, requesterName, utils::MasterSlave::_rank);

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
  Event e2("m2n.buildCommunicationMap", precice::syncMode);
  std::map<int, std::vector<int>> communicationMap = m2n::buildCommunicationMap(
    vertexDistribution, requesterVertexDistribution);
  e2.stop();

// Print `communicationMap'.
#ifdef P2P_LCM_PRINT
  print(communicationMap);
#endif

// Print statistics of `communicationMap'.
#ifdef P2P_LCM_PRINT_STATS
  printCommunicationPartnerCountStats(communicationMap);
  printLocalIndexCountStats(communicationMap);
#endif

#ifdef SuperMUC_WORK
  try {
    auto addressDirectory = _communicationFactory->addressDirectory();

    if (utils::MasterSlave::_masterMode) {
      Event e3("m2n.createDirectories");

      for (int rank = 0; rank < utils::MasterSlave::_size; ++rank) {
        Publisher::createDirectory(addressDirectory + "/" + "." + acceptorName + "-" + _mesh->getName() +
                                   "-" + std::to_string(rank) + ".address");
      }
    }
    utils::Parallel::synchronizeProcesses();
  } catch (...) {
  }
#endif

  Event e4("m2n.createCommunications");
  e4.addData("Connections", communicationMap.size());
  if (communicationMap.empty()) {
    _isConnected = true;
    return;
  }

  _communication = _communicationFactory->newCommunication();

#ifdef SuperMUC_WORK
  Publisher::ScopedPushDirectory spd("." + acceptorName + "-" + _mesh->getName() + "-" +
                                     std::to_string(utils::MasterSlave::_rank) + ".address");
#endif

  // Accept point-to-point connections (as server) between the current acceptor
  // process (in the current participant) with rank `utils::MasterSlave::_rank'
  // and (multiple) requester processes (in the requester participant).
  _communication->acceptConnectionAsServer(acceptorName,
                                           requesterName,
                                           utils::MasterSlave::_rank,
                                           communicationMap.size());

  for (auto const & comMap : communicationMap) {
    int globalRequesterRank = comMap.first;
    auto indices = std::move(communicationMap[globalRequesterRank]);

    _mappings.push_back({globalRequesterRank, std::move(indices), com::PtrRequest(), {}});
  }
  e4.stop();
  _isConnected = true;
}

void PointToPointCommunication::requestConnection(std::string const &acceptorName,
                                                  std::string const &requesterName)
{
  TRACE(acceptorName, requesterName);
  CHECK(not isConnected(), "Already connected!");
  CHECK(utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode,
        "You can only use a point-to-point communication between two participants which both use a master. "
        << "Please use distribution-type gather-scatter instead.");

  mesh::Mesh::VertexDistribution &vertexDistribution = _mesh->getVertexDistribution();
  mesh::Mesh::VertexDistribution  acceptorVertexDistribution;

  if (utils::MasterSlave::_masterMode) {
    Event e0("m2n.exchangeVertexDistribution");
    // Establish connection between participants' master processes.
    auto c = _communicationFactory->newCommunication();
    c->requestConnection(acceptorName, requesterName, 0, 1);

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
  Event e2("m2n.buildCommunicationMap", precice::syncMode);
  std::map<int, std::vector<int>> communicationMap = m2n::buildCommunicationMap(
    vertexDistribution, acceptorVertexDistribution);
  e2.stop();

// Print `communicationMap'.
#ifdef P2P_LCM_PRINT
  print(communicationMap);
#endif

// Print statistics of `communicationMap'.
#ifdef P2P_LCM_PRINT_STATS
  printCommunicationPartnerCountStats(communicationMap);
  printLocalIndexCountStats(communicationMap);
#endif

#ifdef SuperMUC_WORK
  try {
    auto addressDirectory = _communicationFactory->addressDirectory();

    utils::Parallel::synchronizeProcesses();
  } catch (...) {
  }
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

  _communication = _communicationFactory->newCommunication();
  // Request point-to-point connections (as client) between the current
  // requester process (in the current participant) and (multiple) acceptor
  // processes (in the acceptor participant) to ranks `accceptingRanks'
  // according to `communicationMap`.
  _communication->requestConnectionAsClient(acceptorName, requesterName,
                                            acceptingRanks, utils::MasterSlave::_rank);

  for (auto &i : communicationMap) {
    auto globalAcceptorRank = i.first;
    auto indices            = std::move(i.second);

#ifdef SuperMUC_WORK
    Publisher::ScopedPushDirectory spd("." + acceptorName + "-" + _mesh->getName() + "-" +
                                       std::to_string(globalAcceptorRank) + ".address");
#endif

    _mappings.push_back({globalAcceptorRank, std::move(indices), com::PtrRequest(), {}});
  }
  e4.stop();
  _isConnected = true;
}

void PointToPointCommunication::closeConnection()
{
  TRACE();

  if (not isConnected())
    return;

  checkBufferedRequests(true);

  _communication.reset();
  _mappings.clear();
  _isConnected = false;
}

void PointToPointCommunication::send(double *itemsToSend,
                                     size_t  size,
                                     int     valueDimension)
{

  if (_mappings.empty()) {
    return;
  }

  for (auto &mapping : _mappings) {
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

void PointToPointCommunication::checkBufferedRequests(bool blocking)
{
  TRACE(bufferedRequests.size());
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
