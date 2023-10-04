#include <algorithm>
#include <boost/container/flat_map.hpp>
#include <boost/io/ios_state.hpp>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <thread>
#include <utility>
#include <vector>

#include "PointToPointCommunication.hpp"
#include "com/Communication.hpp"
#include "com/CommunicationFactory.hpp"
#include "com/Extra.hpp"
#include "com/Request.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "profiling/Event.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

using precice::profiling::Event;

namespace precice::m2n {

void send(mesh::Mesh::VertexDistribution const &m,
          int                                   rankReceiver,
          const com::PtrCommunication &         communication)
{
  communication->send(static_cast<int>(m.size()), rankReceiver);

  for (auto const &i : m) {
    auto const &rank    = i.first;
    auto const &indices = i.second;
    communication->send(rank, rankReceiver);
    communication->sendRange(indices, rankReceiver);
  }
}

void receive(mesh::Mesh::VertexDistribution &m,
             int                             rankSender,
             const com::PtrCommunication &   communication)
{
  using precice::com::AsVectorTag;
  m.clear();
  int size = 0;
  communication->receive(size, rankSender);

  while (size--) {
    Rank rank = -1;
    communication->receive(rank, rankSender);
    m[rank] = communication->receiveRange(rankSender, AsVectorTag<int>{});
  }
}

void broadcastSend(mesh::Mesh::VertexDistribution const &m,
                   const com::PtrCommunication &         communication = utils::IntraComm::getCommunication())
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
                      const com::PtrCommunication &   communication = utils::IntraComm::getCommunication())
{
  m.clear();
  int size = 0;
  communication->broadcast(size, rankBroadcaster);

  while (size--) {
    Rank rank = -1;
    communication->broadcast(rank, rankBroadcaster);
    communication->broadcast(m[rank], rankBroadcaster);
  }
}

void broadcast(mesh::Mesh::VertexDistribution &m)
{
  if (utils::IntraComm::isPrimary()) {
    // Broadcast (send) vertex distributions.
    m2n::broadcastSend(m);
  } else if (utils::IntraComm::isSecondary()) {
    // Broadcast (receive) vertex distributions.
    m2n::broadcastReceive(m, 0);
  }
}

void print(std::map<int, std::vector<int>> const &m)
{
  std::ostringstream oss;

  oss << "rank: " << utils::IntraComm::getRank() << "\n";

  for (auto &i : m) {
    for (auto &j : i.second) {
      oss << i.first << ":" << j << '\n'; // prints rank:index
    }
  }

  if (utils::IntraComm::isSecondary()) {
    utils::IntraComm::getCommunication()->send(oss.str(), 0);
  } else {

    std::string s;

    for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
      utils::IntraComm::getCommunication()->receive(s, rank);

      oss << s;
    }

    std::cout << oss.str();
  }
}

void printCommunicationPartnerCountStats(std::map<int, std::vector<int>> const &m)
{
  int size = m.size();

  if (utils::IntraComm::isPrimary()) {
    size_t count   = 0;
    size_t maximum = std::numeric_limits<size_t>::min();
    size_t minimum = std::numeric_limits<size_t>::max();
    size_t total   = size;

    if (size) {
      maximum = std::max(maximum, static_cast<size_t>(size));
      minimum = std::min(minimum, static_cast<size_t>(size));
      count++;
    }

    for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
      utils::IntraComm::getCommunication()->receive(size, rank);

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

    boost::io::ios_all_saver ias{std::cout};
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
    PRECICE_ASSERT(utils::IntraComm::isSecondary());
    utils::IntraComm::getCommunication()->send(size, 0);
  }
}

void printLocalIndexCountStats(std::map<int, std::vector<int>> const &m)
{
  int size = 0;

  for (auto &i : m) {
    size += i.second.size();
  }

  if (utils::IntraComm::isPrimary()) {
    size_t count   = 0;
    size_t maximum = std::numeric_limits<size_t>::min();
    size_t minimum = std::numeric_limits<size_t>::max();
    size_t total   = size;

    if (size) {
      maximum = std::max(maximum, static_cast<size_t>(size));
      minimum = std::min(minimum, static_cast<size_t>(size));

      count++;
    }

    for (Rank rank : utils::IntraComm::allSecondaryRanks()) {
      utils::IntraComm::getCommunication()->receive(size, rank);

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

    boost::io::ios_all_saver ias{std::cout};
    std::cout << std::fixed << std::setprecision(3) //
              << "Number of LVDIs per Interface Process:"
              << "\n"
              << "  Total:   " << total << '\n'
              << "  Maximum: " << maximum << '\n'
              << "  Minimum: " << minimum << '\n'
              << "  Average: " << average << '\n'
              << "Number of Interface Processes: " << count << '\n'
              << '\n';
  } else {
    PRECICE_ASSERT(utils::IntraComm::isSecondary());

    utils::IntraComm::getCommunication()->send(size, 0);
  }
}

namespace {
/**
   * @brief This function is by by and large the same as std::set_intersection().
   * The only difference is that we don't return the intersection set itself, but
   * we return the indices of elements in \p InputIt1, which appear in both sets
   * ( \p InputIt1 and \p InputIt2 )
   * The implementation was taken from
   * https://en.cppreference.com/w/cpp/algorithm/set_intersection#Version_1 with the
   * only difference that we compute and store std::distance() in the output iterator.
   * Similar to the std function, this function operates on sorted ranges.
   *
   * @param ref1 The reference iterator, to which we compute the distance/indices.
   * @param first1 The begin of the first range we want to compute the intersection with
   * @param last1 The end of the first range we want to compute the intersection with
   * @param first1 The begin of the second range we want to compute the intersection with
   * @param last1 The end of the second range we want to compute the intersection with
   * @param d_first Beginning of the output range
   * @return OutputIt Iterator past the end of the constructed range.
   */
template <class InputIt1, class InputIt2, class OutputIt>
OutputIt set_intersection_indices(InputIt1 ref1, InputIt1 first1, InputIt1 last1,
                                  InputIt2 first2, InputIt2 last2, OutputIt d_first)
{
  while (first1 != last1 && first2 != last2) {
    if (*first1 < *first2)
      ++first1;
    else {
      if (!(*first2 < *first1))
        *d_first++ = std::distance(ref1, first1++); // *first1 and *first2 are equivalent.
      ++first2;
    }
  }
  return d_first;
}
} // namespace

/** builds the communication map for a local distribution given the global distribution.
 *
 *
 * @param[in] thisVertexDistribution the local vertex distribution
 * @param[in] otherVertexDistribution the total vertex distribution
 * @param[in] thisRank the rank to build the map for
 *
 * @returns the resulting communication map for rank thisRank
 *
 * The worst case complexity of the function is:
 * \f$ \mathcal{O}(p n \log(n) + p 2 (2 n)) \f$
 *
 * which is composed of the initial std::sort for each vector and the subsequent
 * computation of the intersection.
 *
 * * n is the number of data indices for each vector in `otherVertexDistribution'
 * * p number of ranks
 * * Note that n becomes smaller, if we have more ranks.
 *
 * However, in case our indices are already sorted (which is typically the case),
 * the complexity boils down to linear complexity
 * \f$ \mathcal{O}(p n) \f$
 */
std::map<int, std::vector<int>> buildCommunicationMap(
    // `thisVertexDistribution' is input vertex distribution from this participant.
    mesh::Mesh::VertexDistribution &thisVertexDistribution,
    // `otherVertexDistribution' is input vertex distribution from other participant.
    mesh::Mesh::VertexDistribution &otherVertexDistribution,
    int                             thisRank = utils::IntraComm::getRank())
{
  auto iterator = thisVertexDistribution.find(thisRank);
  if (iterator == thisVertexDistribution.end())
    return {};

  std::map<int, std::vector<int>> communicationMap;

  // take advantage that these data structures are in most cases sorted by construction,
  // i.e., we perform here mostly a safety check and don't perform an actual sorting
  if (!std::is_sorted(iterator->second.begin(), iterator->second.end())) {
    std::sort(iterator->second.begin(), iterator->second.end());
  }

  // now we iterate over all other vertex distributions to compute the intersection
  for (auto &other : otherVertexDistribution) {
    // first a safety check, that we are actually sorted, similar to above
    if (!std::is_sorted(other.second.begin(), other.second.end())) {
      std::sort(other.second.begin(), other.second.end());
    }

    // before starting to compute an actual intersection, we first check if elements can
    // possibly be in both data sets by comparing upper and lower index bounds of both
    // data sets. For typical partitioning schemes, each rank only exchanges data with
    // a few neighbors such that this check already filters out a significant amount of
    // computations
    if (iterator->second.empty() || other.second.empty() || (other.second.back() < iterator->second.at(0)) || (other.second.at(0) > iterator->second.back())) {
      // in this case there is nothing to be done
      continue;
    } else {
      // we have an intersection, let's compute it
      std::vector<int> inters;
      // the actual worker function, which gives us the indices of intersecting elements
      // have a look at the documentation of the function for more details
      precice::m2n::set_intersection_indices(iterator->second.begin(), iterator->second.begin(), iterator->second.end(),
                                             other.second.begin(), other.second.end(),
                                             std::back_inserter(inters));
      // we have the results, now commit it into the final map
      if (!inters.empty())
        communicationMap.insert({other.first, std::move(inters)});
    }
  }
  return communicationMap;
}

PointToPointCommunication::PointToPointCommunication(
    com::PtrCommunicationFactory communicationFactory,
    mesh::PtrMesh                mesh)
    : DistributedCommunication(std::move(mesh)),
      _communicationFactory(std::move(communicationFactory))
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

  mesh::Mesh::VertexDistribution vertexDistribution = _mesh->getVertexDistribution();
  mesh::Mesh::VertexDistribution requesterVertexDistribution;

  if (not utils::IntraComm::isSecondary()) {
    PRECICE_DEBUG("Exchange vertex distribution between both primary ranks");
    Event e0("m2n.exchangeVertexDistribution");
    // Establish connection between participants' primary processes.
    auto c = _communicationFactory->newCommunication();

    c->acceptConnection(acceptorName, requesterName, "TMP-PRIMARYCOM-" + _mesh->getName(), utils::IntraComm::getRank());

    // Exchange vertex distributions.
    m2n::send(vertexDistribution, 0, c);
    m2n::receive(requesterVertexDistribution, 0, c);
  }

  PRECICE_DEBUG("Broadcast vertex distributions");
  Event e1("m2n.broadcastVertexDistributions", profiling::Synchronize);
  m2n::broadcast(vertexDistribution);
  if (utils::IntraComm::isSecondary()) {
    _mesh->setVertexDistribution(vertexDistribution);
  }
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
  Event                           e2("m2n.buildCommunicationMap", profiling::Synchronize);
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
  // process (in the current participant) with rank `utils::IntraComm::getRank()'
  // and (multiple) requester processes (in the requester participant).
  _communication->acceptConnectionAsServer(acceptorName,
                                           requesterName,
                                           _mesh->getName(),
                                           utils::IntraComm::getRank(),
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
      utils::IntraComm::getRank(),
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

  mesh::Mesh::VertexDistribution vertexDistribution = _mesh->getVertexDistribution();
  mesh::Mesh::VertexDistribution acceptorVertexDistribution;

  if (not utils::IntraComm::isSecondary()) {
    PRECICE_DEBUG("Exchange vertex distribution between both primary ranks");
    Event e0("m2n.exchangeVertexDistribution");
    // Establish connection between participants' primary processes.
    auto c = _communicationFactory->newCommunication();
    c->requestConnection(acceptorName, requesterName,
                         "TMP-PRIMARYCOM-" + _mesh->getName(),
                         0, 1);

    // Exchange vertex distributions.
    m2n::receive(acceptorVertexDistribution, 0, c);
    m2n::send(vertexDistribution, 0, c);
  }

  PRECICE_DEBUG("Broadcast vertex distributions");
  Event e1("m2n.broadcastVertexDistributions", profiling::Synchronize);
  m2n::broadcast(vertexDistribution);
  if (utils::IntraComm::isSecondary()) {
    _mesh->setVertexDistribution(vertexDistribution);
  }
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
  Event                           e2("m2n.buildCommunicationMap", profiling::Synchronize);
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
                                            acceptingRanks, utils::IntraComm::getRank());

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
                                            acceptingRanks, utils::IntraComm::getRank());

  for (auto &connectedRank : localConnectedRanks) {
    _connectionDataVector.push_back({connectedRank, com::PtrRequest()});
  }
  _isConnected = true;
}

void PointToPointCommunication::completeSecondaryRanksConnection()
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
  _connectionDataVector.clear();
  _isConnected = false;
}

void PointToPointCommunication::send(precice::span<double const> itemsToSend, int valueDimension)
{

  if (_mappings.empty() || itemsToSend.empty()) {
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
    auto request = _communication->aSend(span<const double>{*buffer}, mapping.remoteRank);
    bufferedRequests.emplace_back(request, buffer);
  }
  checkBufferedRequests(false);
}

void PointToPointCommunication::receive(precice::span<double> itemsToReceive, int valueDimension)
{
  if (_mappings.empty() || itemsToReceive.empty()) {
    return;
  }

  std::fill(itemsToReceive.begin(), itemsToReceive.end(), 0.0);

  for (auto &mapping : _mappings) {
    mapping.recvBuffer.resize(mapping.indices.size() * valueDimension);
    mapping.request = _communication->aReceive(span<double>{mapping.recvBuffer}, mapping.remoteRank);
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

void PointToPointCommunication::broadcastSend(int itemToSend)
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
    com::sendMesh(*_communication, connectionData.remoteRank, *_mesh);
  }
}

void PointToPointCommunication::broadcastReceiveAllMesh()
{
  for (auto &connectionData : _connectionDataVector) {
    com::receiveMesh(*_communication, connectionData.remoteRank, *_mesh);
  }
}

void PointToPointCommunication::scatterAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  for (auto &connectionData : _connectionDataVector) {
    _communication->sendRange(localCommunicationMap[connectionData.remoteRank], connectionData.remoteRank);
  }
}

void PointToPointCommunication::gatherAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  using precice::com::AsVectorTag;
  for (auto &connectionData : _connectionDataVector) {
    localCommunicationMap[connectionData.remoteRank] = _communication->receiveRange(connectionData.remoteRank, AsVectorTag<int>{});
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

} // namespace precice::m2n
