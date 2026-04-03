#include "com/SocketIntraComm.hpp"
#include "com/SocketCommunication.hpp"
#include "logging/LogMacros.hpp"
#include "precice/impl/Types.hpp"
#include "utils/assertion.hpp"
#include "utils/networking.hpp"
#include "utils/span_tools.hpp"

#include <algorithm>
#include <memory>
#include <ostream>
#include <vector>

namespace precice::com {

SocketIntraComm::SocketIntraComm(unsigned short portNumber,
                                 bool           reuseAddress,
                                 std::string    networkName,
                                 std::string    addressDirectory)
    : _socketComm(std::make_unique<SocketCommunication>(
          portNumber, reuseAddress,
          networkName.empty() ? utils::networking::loopbackInterfaceName() : std::move(networkName),
          std::move(addressDirectory)))
{
}

SocketIntraComm::~SocketIntraComm()
{
  PRECICE_TRACE(_isConnected);
  closeConnection();
}

bool SocketIntraComm::isConnected()
{
  return _isConnected;
}

size_t SocketIntraComm::getRemoteCommunicatorSize()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(isConnected());
  return _remoteCommunicatorSize;
}

void SocketIntraComm::connectIntraComm(std::string const &participantName,
                                       std::string const &tag,
                                       int                rank,
                                       int                size)
{
  PRECICE_TRACE(participantName, rank, size);
  if (size == 1)
    return;

  // Delegate connection establishment to the internal SocketCommunication
  // using the same pattern as Communication::connectIntraComm
  std::string primaryName   = participantName + "Primary";
  std::string secondaryName = participantName + "Secondary";

  int secondaryRanksSize = size - RANK_OFFSET;
  if (rank == 0) {
    PRECICE_INFO("Connecting Primary rank to {} Secondary ranks via sockets", secondaryRanksSize);
    _socketComm->prepareEstablishment(primaryName, secondaryName);
    _socketComm->acceptConnection(primaryName, secondaryName, tag, rank, RANK_OFFSET);
    _socketComm->cleanupEstablishment(primaryName, secondaryName);
  } else {
    int secondaryRank = rank - RANK_OFFSET;
    PRECICE_INFO("Connecting Secondary rank #{} to Primary rank via sockets", secondaryRank);
    _socketComm->requestConnection(primaryName, secondaryName, tag, secondaryRank, secondaryRanksSize);
  }

  _remoteCommunicatorSize = _socketComm->getRemoteCommunicatorSize();
  _isConnected            = true;
}

void SocketIntraComm::closeConnection()
{
  PRECICE_TRACE();
  if (not isConnected())
    return;
  _socketComm->closeConnection();
  _isConnected = false;
}

// ==========================================================================
// Reduction — built from point-to-point socket operations
// ==========================================================================

void SocketIntraComm::reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());

  std::copy(itemsToSend.begin(), itemsToSend.end(), itemsToReceive.begin());

  std::vector<double> received(itemsToReceive.size());
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    auto request = _socketComm->aReceive(received, static_cast<int>(rank) + RANK_OFFSET);
    request->wait();
    for (size_t i = 0; i < itemsToReceive.size(); i++) {
      itemsToReceive[i] += received[i];
    }
  }
}

void SocketIntraComm::reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());
  // Secondary side: primaryRank maps directly to _sockets[0] (no offset needed)
  auto request = _socketComm->aSend(itemsToSend, primaryRank);
  request->wait();
}

void SocketIntraComm::reduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();
  itemToReceive = itemToSend;
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    auto request = _socketComm->aReceive(itemToSend, static_cast<int>(rank) + RANK_OFFSET);
    request->wait();
    itemToReceive += itemToSend;
  }
}

void SocketIntraComm::reduceSum(int itemToSend, int &itemToReceive, Rank primaryRank)
{
  PRECICE_TRACE();
  auto request = _socketComm->aSend(precice::refToSpan<const int>(itemToSend), primaryRank);
  request->wait();
}

void SocketIntraComm::allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());

  reduceSum(itemsToSend, itemsToReceive);

  std::vector<PtrRequest> requests;
  requests.reserve(_remoteCommunicatorSize);
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    requests.push_back(_socketComm->aSend(itemsToReceive, static_cast<int>(rank) + RANK_OFFSET));
  }
  Request::wait(requests);
}

void SocketIntraComm::allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank)
{
  PRECICE_TRACE(itemsToSend.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());

  reduceSum(itemsToSend, itemsToReceive, primaryRank);
  _socketComm->receive(itemsToReceive, primaryRank);
}

void SocketIntraComm::allreduceSum(double itemToSend, double &itemToReceive)
{
  PRECICE_TRACE();
  itemToReceive = itemToSend;
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    auto request = _socketComm->aReceive(itemToSend, static_cast<int>(rank) + RANK_OFFSET);
    request->wait();
    itemToReceive += itemToSend;
  }
  std::vector<PtrRequest> requests(_remoteCommunicatorSize);
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    requests[rank] = _socketComm->aSend(itemToReceive, static_cast<int>(rank) + RANK_OFFSET);
  }
  Request::wait(requests);
}

void SocketIntraComm::allreduceSum(double itemToSend, double &itemToReceive, Rank primaryRank)
{
  PRECICE_TRACE();
  auto request = _socketComm->aSend(itemToSend, primaryRank);
  request->wait();
  _socketComm->receive(itemToReceive, primaryRank);
}

void SocketIntraComm::allreduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();
  itemToReceive = itemToSend;
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    auto request = _socketComm->aReceive(itemToSend, static_cast<int>(rank) + RANK_OFFSET);
    request->wait();
    itemToReceive += itemToSend;
  }
  std::vector<PtrRequest> requests(_remoteCommunicatorSize);
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    requests[rank] = _socketComm->aSend(itemToReceive, static_cast<int>(rank) + RANK_OFFSET);
  }
  Request::wait(requests);
}

void SocketIntraComm::allreduceSum(int itemToSend, int &itemToReceive, Rank primaryRank)
{
  PRECICE_TRACE();
  auto request = _socketComm->aSend(precice::refToSpan<const int>(itemToSend), primaryRank);
  request->wait();
  _socketComm->receive(itemToReceive, primaryRank);
}

// ==========================================================================
// Broadcast — built from point-to-point socket operations
// ==========================================================================

void SocketIntraComm::broadcast(precice::span<const int> itemsToSend)
{
  PRECICE_TRACE(itemsToSend.size());
  std::vector<PtrRequest> requests(_remoteCommunicatorSize);
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    requests[rank] = _socketComm->aSend(itemsToSend, static_cast<int>(rank) + RANK_OFFSET);
  }
  Request::wait(requests);
}

void SocketIntraComm::broadcast(precice::span<int> itemsToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE(itemsToReceive.size());
  _socketComm->receive(itemsToReceive, rankBroadcaster);
}

void SocketIntraComm::broadcast(int itemToSend)
{
  PRECICE_TRACE();
  std::vector<PtrRequest> requests(_remoteCommunicatorSize);
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    requests[rank] = _socketComm->aSend(itemToSend, static_cast<int>(rank) + RANK_OFFSET);
  }
  Request::wait(requests);
}

void SocketIntraComm::broadcast(int &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  _socketComm->receive(itemToReceive, rankBroadcaster);
}

void SocketIntraComm::broadcast(precice::span<const double> itemsToSend)
{
  PRECICE_TRACE(itemsToSend.size());
  std::vector<PtrRequest> requests(_remoteCommunicatorSize);
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    requests[rank] = _socketComm->aSend(itemsToSend, static_cast<int>(rank) + RANK_OFFSET);
  }
  Request::wait(requests);
}

void SocketIntraComm::broadcast(precice::span<double> itemsToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE(itemsToReceive.size());
  _socketComm->receive(itemsToReceive, rankBroadcaster);
}

void SocketIntraComm::broadcast(double itemToSend)
{
  PRECICE_TRACE();
  std::vector<PtrRequest> requests(_remoteCommunicatorSize);
  for (size_t rank = 0; rank < _remoteCommunicatorSize; ++rank) {
    requests[rank] = _socketComm->aSend(itemToSend, static_cast<int>(rank) + RANK_OFFSET);
  }
  Request::wait(requests);
}

void SocketIntraComm::broadcast(double &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  _socketComm->receive(itemToReceive, rankBroadcaster);
}

void SocketIntraComm::broadcast(bool itemToSend)
{
  PRECICE_TRACE();
  int item = itemToSend;
  broadcast(item);
}

void SocketIntraComm::broadcast(bool &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  int item;
  broadcast(item, rankBroadcaster);
  itemToReceive = item;
}

void SocketIntraComm::broadcast(std::vector<int> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(precice::span<const int>{v});
}

void SocketIntraComm::broadcast(std::vector<int> &v, Rank rankBroadcaster)
{
  int size = 0;
  broadcast(size, rankBroadcaster);
  v.clear();
  v.resize(size);
  broadcast(precice::span<int>{v}, rankBroadcaster);
}

void SocketIntraComm::broadcast(std::vector<double> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(precice::span<const double>{v});
}

void SocketIntraComm::broadcast(std::vector<double> &v, Rank rankBroadcaster)
{
  int size = 0;
  broadcast(size, rankBroadcaster);
  v.clear();
  v.resize(size);
  broadcast(precice::span<double>{v}, rankBroadcaster);
}

// ==========================================================================
// Point-to-point — delegated to internal SocketCommunication
// Rank adjustment: IntraCommunication ranks are direct participant ranks
// (0=primary, 1..N-1=secondaries). The internal SocketCommunication uses
// per-rank offset, so we delegate directly.
// ==========================================================================

void SocketIntraComm::send(std::string const &itemToSend, Rank rankReceiver)
{
  _socketComm->send(itemToSend, rankReceiver);
}

void SocketIntraComm::send(precice::span<const int> itemsToSend, Rank rankReceiver)
{
  _socketComm->send(itemsToSend, rankReceiver);
}

PtrRequest SocketIntraComm::aSend(precice::span<const int> itemsToSend, Rank rankReceiver)
{
  return _socketComm->aSend(itemsToSend, rankReceiver);
}

void SocketIntraComm::send(precice::span<const double> itemsToSend, Rank rankReceiver)
{
  _socketComm->send(itemsToSend, rankReceiver);
}

PtrRequest SocketIntraComm::aSend(precice::span<const double> itemsToSend, Rank rankReceiver)
{
  return _socketComm->aSend(itemsToSend, rankReceiver);
}

void SocketIntraComm::send(double itemToSend, Rank rankReceiver)
{
  _socketComm->send(itemToSend, rankReceiver);
}

PtrRequest SocketIntraComm::aSend(const double &itemToSend, Rank rankReceiver)
{
  return _socketComm->aSend(itemToSend, rankReceiver);
}

void SocketIntraComm::send(int itemToSend, Rank rankReceiver)
{
  _socketComm->send(itemToSend, rankReceiver);
}

PtrRequest SocketIntraComm::aSend(const int &itemToSend, Rank rankReceiver)
{
  return _socketComm->aSend(itemToSend, rankReceiver);
}

void SocketIntraComm::send(bool itemToSend, Rank rankReceiver)
{
  _socketComm->send(itemToSend, rankReceiver);
}

PtrRequest SocketIntraComm::aSend(const bool &itemToSend, Rank rankReceiver)
{
  return _socketComm->aSend(itemToSend, rankReceiver);
}

void SocketIntraComm::receive(std::string &itemToReceive, Rank rankSender)
{
  _socketComm->receive(itemToReceive, rankSender);
}

void SocketIntraComm::receive(precice::span<int> itemsToReceive, Rank rankSender)
{
  _socketComm->receive(itemsToReceive, rankSender);
}

void SocketIntraComm::receive(precice::span<double> itemsToReceive, Rank rankSender)
{
  _socketComm->receive(itemsToReceive, rankSender);
}

PtrRequest SocketIntraComm::aReceive(precice::span<double> itemsToReceive, int rankSender)
{
  return _socketComm->aReceive(itemsToReceive, rankSender);
}

void SocketIntraComm::receive(double &itemToReceive, Rank rankSender)
{
  _socketComm->receive(itemToReceive, rankSender);
}

PtrRequest SocketIntraComm::aReceive(double &itemToReceive, Rank rankSender)
{
  return _socketComm->aReceive(itemToReceive, rankSender);
}

void SocketIntraComm::receive(int &itemToReceive, Rank rankSender)
{
  _socketComm->receive(itemToReceive, rankSender);
}

PtrRequest SocketIntraComm::aReceive(int &itemToReceive, Rank rankSender)
{
  return _socketComm->aReceive(itemToReceive, rankSender);
}

void SocketIntraComm::receive(bool &itemToReceive, Rank rankSender)
{
  _socketComm->receive(itemToReceive, rankSender);
}

PtrRequest SocketIntraComm::aReceive(bool &itemToReceive, Rank rankSender)
{
  return _socketComm->aReceive(itemToReceive, rankSender);
}

} // namespace precice::com
