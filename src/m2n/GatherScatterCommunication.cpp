#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <ostream>
#include <utility>

#include "GatherScatterCommunication.hpp"
#include "com/Communication.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "mesh/Mesh.hpp"
#include "precice/impl/Types.hpp"
#include "utils/IntraComm.hpp"
#include "utils/algorithm.hpp"
#include "utils/assertion.hpp"

namespace precice::m2n {
GatherScatterCommunication::GatherScatterCommunication(
    com::PtrCommunication com,
    mesh::PtrMesh         mesh)
    : DistributedCommunication(std::move(mesh)),
      _com(std::move(com)),
      _isConnected(false)
{
}

GatherScatterCommunication::~GatherScatterCommunication()
{
  if (isConnected()) {
    closeConnection();
  }
}

bool GatherScatterCommunication::isConnected() const
{
  return _isConnected;
}

void GatherScatterCommunication::acceptConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(utils::IntraComm::isSecondary() || _com->isConnected());
  _isConnected = true;
}

void GatherScatterCommunication::requestConnection(
    const std::string &acceptorName,
    const std::string &requesterName)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(utils::IntraComm::isSecondary() || _com->isConnected());
  _isConnected = true;
}

void GatherScatterCommunication::closeConnection()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(utils::IntraComm::isSecondary() || not _com->isConnected());
  _isConnected = false;
}

namespace {
template <class Indices, class Src, class Dst, class Size>
void add_to_indirect_blocks(
    const Src     &src,
    const Indices &indices,
    Size           blockSize,
    Dst           &dst)
{
  for (size_t i = 0; i < indices.size(); ++i) {
    auto srcfirst = blockSize * i;
    auto dstfirst = blockSize * indices[i];
    PRECICE_ASSERT(srcfirst >= 0 && static_cast<size_t>(srcfirst + blockSize) <= src.size(), srcfirst, blockSize, src.size());
    PRECICE_ASSERT(dstfirst >= 0 && static_cast<size_t>(dstfirst + blockSize) <= dst.size(), dstfirst, blockSize, dst.size());
    utils::add_n(&src[srcfirst], blockSize, &dst[dstfirst]);
  }
}

template <class Indices, class Src, class Dst, class Size>
void copy_from_indirect_blocks(
    const Src     &src,
    const Indices &indices,
    Size           blockSize,
    Dst           &dst)
{
  for (size_t i = 0; i < indices.size(); ++i) {
    auto srcfirst = blockSize * indices[i];
    auto dstfirst = blockSize * i;
    PRECICE_ASSERT(srcfirst >= 0 && static_cast<size_t>(srcfirst + blockSize) <= src.size(), srcfirst, blockSize, src.size());
    PRECICE_ASSERT(dstfirst >= 0 && static_cast<size_t>(dstfirst + blockSize) <= dst.size(), dstfirst, blockSize, dst.size());
    std::copy_n(&src[srcfirst], blockSize, &dst[dstfirst]);
  }
}

} // namespace

void GatherScatterCommunication::send(precice::span<double const> itemsToSend, int valueDimension)
{
  PRECICE_TRACE(itemsToSend.size());

  // Gather data on secondary ranks
  if (utils::IntraComm::isSecondary()) { // Secondary rank
    if (!itemsToSend.empty()) {
      PRECICE_DEBUG("Providing {} elements to gather step", itemsToSend.size());
      utils::IntraComm::getCommunication()->send(itemsToSend, 0);
    }
    return;
  }

  // Primary rank or coupling mode
  PRECICE_ASSERT(utils::IntraComm::getRank() == 0);
  const auto &vertexDistribution = _mesh->getVertexDistribution();
  const int   globalSize         = _mesh->getGlobalNumberOfVertices() * valueDimension;
  PRECICE_DEBUG("Gathering data on primary ({} elements)", globalSize);
  std::vector<double> globalItemsToSend(globalSize);

  // Directly copy primary rank data
  PRECICE_ASSERT(vertexDistribution.count(0) > 0);
  const auto &primaryDistribution = vertexDistribution.at(0);
  add_to_indirect_blocks(itemsToSend, primaryDistribution, valueDimension, globalItemsToSend);
  PRECICE_DEBUG("Directly gathered {} entries from primary", primaryDistribution.size() * valueDimension);

  // Gather data from secondary ranks
  for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
    PRECICE_ASSERT(utils::IntraComm::getCommunication() != nullptr);
    PRECICE_ASSERT(utils::IntraComm::getCommunication()->isConnected());

    auto iter = vertexDistribution.find(secondaryRank);
    if (iter == vertexDistribution.end()) {
      continue;
    }
    const auto &secondaryDistribution = iter->second;

    int secondaryRankSize = secondaryDistribution.size() * valueDimension;
    PRECICE_DEBUG("Gathering {} entries from secondary rank {}", secondaryRankSize, secondaryRank);
    if (secondaryDistribution.empty()) {
      continue;
    }
    std::vector<double> secondaryRankValues(secondaryRankSize);
    utils::IntraComm::getCommunication()->receive(span<double>{secondaryRankValues}, secondaryRank);
    add_to_indirect_blocks(secondaryRankValues, secondaryDistribution, valueDimension, globalItemsToSend);
  }

  // Send data to other primary
  PRECICE_DEBUG("Sending gathered data to other participant");
  _com->sendRange(globalItemsToSend, 0);
}

void GatherScatterCommunication::receive(precice::span<double> itemsToReceive, int valueDimension)
{
  PRECICE_TRACE(itemsToReceive.size());

  // Secondary ranks receive scattered data
  if (utils::IntraComm::isSecondary()) { // Secondary rank
    if (!itemsToReceive.empty()) {
      auto received = utils::IntraComm::getCommunication()->receiveRange(0, com::asVector<double>);
      PRECICE_ASSERT(!received.empty());
      PRECICE_DEBUG("Received scattered data starting with {}", received[0]);
      std::copy(received.begin(), received.end(), itemsToReceive.begin());
    }
    return;
  }

  // Primary rank receives and scatters the data
  PRECICE_ASSERT(not utils::IntraComm::isSecondary());

  const int globalSize = _mesh->getGlobalNumberOfVertices() * valueDimension;
  PRECICE_DEBUG("Receiving {} elements from other participant to scatter", globalSize);

  auto globalItemsToReceive = _com->receiveRange(0, com::asVector<double>);
  PRECICE_ASSERT(globalItemsToReceive.size() == static_cast<std::size_t>(globalSize));

  const auto &vertexDistribution = _mesh->getVertexDistribution();

  // Directly copy primary rank data
  PRECICE_ASSERT(vertexDistribution.count(0) > 0);
  const auto &primaryDistribution = vertexDistribution.at(0);
  copy_from_indirect_blocks(globalItemsToReceive, primaryDistribution, valueDimension, itemsToReceive);

  PRECICE_DEBUG("Directly extracted {} data entries for primary", primaryDistribution.size() * valueDimension);

  // Extract and scatter data to secondary ranks
  for (Rank secondaryRank : utils::IntraComm::allSecondaryRanks()) {
    PRECICE_ASSERT(utils::IntraComm::getCommunication() != nullptr);
    PRECICE_ASSERT(utils::IntraComm::getCommunication()->isConnected());

    auto iter = vertexDistribution.find(secondaryRank);
    if (iter == vertexDistribution.end()) {
      continue;
    }
    const auto &secondaryDistribution = iter->second;

    int secondarySize = secondaryDistribution.size() * valueDimension;
    PRECICE_DEBUG("Scattering {} entries to secondary {}", secondarySize, secondarySize);
    if (secondaryDistribution.empty()) {
      continue;
    }

    std::vector<double> secondaryRankValues(secondarySize);
    copy_from_indirect_blocks(globalItemsToReceive, secondaryDistribution, valueDimension, secondaryRankValues);
    PRECICE_DEBUG("Scattering data starting with {} to rank {}", secondaryRankValues[0], secondaryRank);
    utils::IntraComm::getCommunication()->sendRange(secondaryRankValues, secondaryRank);
  }
}

void GatherScatterCommunication::acceptPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::requestPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::broadcastSend(int itemToSend)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::broadcastReceiveAll(std::vector<int> &itemToReceive)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::broadcastSendMesh()
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::broadcastReceiveAllMesh()
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::scatterAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::gatherAllCommunicationMap(CommunicationMap &localCommunicationMap)
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

void GatherScatterCommunication::completeSecondaryRanksConnection()
{
  PRECICE_ASSERT(false, "Not available for GatherScatterCommunication.");
}

} // namespace precice::m2n
