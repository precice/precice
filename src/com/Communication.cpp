#include <algorithm>
#include <memory>
#include <ostream>
#include <vector>

#include "Communication.hpp"
#include "Request.hpp"
#include "logging/LogMacros.hpp"
#include "precice/impl/Types.hpp"
#include "utils/assertion.hpp"

namespace precice::com {

void Communication::connectIntraComm(std::string const &participantName,
                                     std::string const &tag,
                                     int                rank,
                                     int                size)
{
  if (size == 1)
    return;

  std::string primaryName   = participantName + "Primary";
  std::string secondaryName = participantName + "Secondary";

  constexpr Rank rankOffset         = 1;
  int            secondaryRanksSize = size - rankOffset;
  if (rank == 0) {
    PRECICE_INFO("Connecting Primary rank to {} Secondary ranks", secondaryRanksSize);
    prepareEstablishment(primaryName, secondaryName);
    acceptConnection(primaryName, secondaryName, tag, rank, rankOffset);
    cleanupEstablishment(primaryName, secondaryName);
  } else {
    int secondaryRank = rank - rankOffset;
    PRECICE_INFO("Connecting Secondary rank #{} to Primary rank", secondaryRank);
    requestConnection(primaryName, secondaryName, tag, secondaryRank, secondaryRanksSize);
  }
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive)
{
  PRECICE_TRACE(itemsToSend.size(), itemsToReceive.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());

  std::copy(itemsToSend.begin(), itemsToSend.end(), itemsToReceive.begin());

  std::vector<double> received(itemsToReceive.size());
  // receive local results from secondary ranks
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request = aReceive(received, rank + _rankOffset);
    request->wait();
    for (size_t i = 0; i < itemsToReceive.size(); i++) {
      itemsToReceive[i] += received[i];
    }
  }
}

void Communication::reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank)
{
  PRECICE_TRACE(itemsToSend.size(), itemsToReceive.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());

  auto request = aSend(itemsToSend, primaryRank);
  request->wait();
}

void Communication::reduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();

  itemToReceive = itemToSend;

  // receive local results from secondary ranks
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }
}

void Communication::reduceSum(int itemToSend, int &itemToReceive, Rank primaryRank)
{
  PRECICE_TRACE();

  auto request = aSend(itemToSend, primaryRank);
  request->wait();
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive)
{
  PRECICE_TRACE(itemsToSend.size(), itemsToReceive.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());

  reduceSum(itemsToSend, itemsToReceive);

  // send reduced result to all secondary ranks
  std::vector<PtrRequest> requests;
  requests.reserve(getRemoteCommunicatorSize());
  for (Rank rank : remoteCommunicatorRanks()) {
    requests.push_back(aSend(itemsToReceive, rank + _rankOffset));
  }
  Request::wait(requests);
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank)
{
  PRECICE_TRACE(itemsToSend.size(), itemsToReceive.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());

  reduceSum(itemsToSend, itemsToReceive, primaryRank);
  // receive reduced data from primary rank
  receive(itemsToReceive, primaryRank + _rankOffset);
}

void Communication::allreduceSum(double itemToSend, double &itemToReceive)
{
  PRECICE_TRACE();

  itemToReceive = itemToSend;

  // receive local results from secondary ranks
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }

  // send reduced result to all secondary ranks
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (Rank rank : remoteCommunicatorRanks()) {
    requests[rank] = aSend(itemToReceive, rank + _rankOffset);
  }
  Request::wait(requests);
}

void Communication::allreduceSum(double itemToSend, double &itemsToReceive, Rank primaryRank)
{
  PRECICE_TRACE();

  auto request = aSend(itemToSend, primaryRank);
  request->wait();
  // receive reduced data from primary rank
  receive(itemsToReceive, primaryRank + _rankOffset);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();

  itemToReceive = itemToSend;

  // receive local results from secondary ranks
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }

  // send reduced result to all secondary ranks
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (Rank rank : remoteCommunicatorRanks()) {
    requests[rank] = aSend(itemToReceive, rank + _rankOffset);
  }
  Request::wait(requests);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive, Rank primaryRank)
{
  PRECICE_TRACE();

  auto request = aSend(itemToSend, primaryRank);
  request->wait();
  // receive reduced data from primary rank
  receive(itemToReceive, primaryRank + _rankOffset);
}

void Communication::broadcast(precice::span<const int> itemsToSend)
{
  PRECICE_TRACE(itemsToSend.size());

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (Rank rank : remoteCommunicatorRanks()) {
    requests[rank] = aSend(itemsToSend, rank + _rankOffset);
  }

  Request::wait(requests);
}

void Communication::broadcast(precice::span<int> itemsToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE(itemsToReceive.size());

  receive(itemsToReceive, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(int itemToSend)
{
  PRECICE_TRACE();

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (Rank rank : remoteCommunicatorRanks()) {
    requests[rank] = aSend(itemToSend, rank + _rankOffset);
  }

  Request::wait(requests);
}

void Communication::broadcast(int &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  receive(itemToReceive, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(precice::span<const double> itemsToSend)
{
  PRECICE_TRACE(itemsToSend.size());

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (Rank rank : remoteCommunicatorRanks()) {
    requests[rank] = aSend(itemsToSend, rank + _rankOffset);
  }

  Request::wait(requests);
}

void Communication::broadcast(precice::span<double> itemsToReceive,
                              int                   rankBroadcaster)
{
  PRECICE_TRACE(itemsToReceive.size());
  receive(itemsToReceive, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(double itemToSend)
{
  PRECICE_TRACE();

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (Rank rank : remoteCommunicatorRanks()) {
    requests[rank] = aSend(itemToSend, rank + _rankOffset);
  }

  Request::wait(requests);
}

void Communication::broadcast(double &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  receive(itemToReceive, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(bool itemToSend)
{
  PRECICE_TRACE();
  int item = itemToSend;
  broadcast(item);
}

void Communication::broadcast(bool &itemToReceive, Rank rankBroadcaster)
{
  PRECICE_TRACE();
  int item;
  broadcast(item, rankBroadcaster);
  itemToReceive = item;
}

void Communication::broadcast(std::vector<int> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(precice::span<const int>{v});
}

void Communication::broadcast(std::vector<int> &v, Rank rankBroadcaster)
{
  int size = 0;
  broadcast(size, rankBroadcaster);

  v.clear();
  v.resize(size);
  broadcast(precice::span<int>{v}, rankBroadcaster);
}

void Communication::broadcast(std::vector<double> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(precice::span<const double>{v});
}

void Communication::broadcast(std::vector<double> &v, Rank rankBroadcaster)
{
  int size = 0;
  broadcast(size, rankBroadcaster);

  v.clear();
  v.resize(size);
  broadcast(precice::span<double>{v}, rankBroadcaster);
}

void Communication::sendRange(precice::span<const double> itemsToSend, Rank rankReceiver)
{
  int size = itemsToSend.size();
  send(size, rankReceiver);
  if (size > 0) {
    send(itemsToSend, rankReceiver);
  }
}

void Communication::sendRange(precice::span<const int> itemsToSend, Rank rankReceiver)
{
  int size = itemsToSend.size();
  send(size, rankReceiver);
  if (size > 0) {
    send(itemsToSend, rankReceiver);
  }
}

std::vector<int> Communication::receiveRange(Rank rankSender, AsVectorTag<int>)
{
  int size{-1};
  receive(size, rankSender);
  PRECICE_ASSERT(size >= 0);
  std::vector<int> result;
  if (size > 0) {
    result.resize(size);
    receive(result, rankSender);
  }
  return result;
}

std::vector<double> Communication::receiveRange(Rank rankSender, AsVectorTag<double>)
{
  int size{-1};
  receive(size, rankSender);
  PRECICE_ASSERT(size >= 0);
  std::vector<double> result;
  if (size > 0) {
    result.resize(size);
    receive(result, rankSender);
  }
  return result;
}

int Communication::adjustRank(Rank rank) const
{
  return rank - _rankOffset;
}

void connectCircularComm(
    std::string const  &participantName,
    std::string const  &tag,
    int                 rank,
    int                 size,
    com::Communication &left,
    com::Communication &right)
{
  PRECICE_ASSERT(!left.isConnected());
  PRECICE_ASSERT(!right.isConnected());
  PRECICE_ASSERT(rank >= 0 && rank < size && size > 0);

  if (size == 1) {
    return;
  }

  const int prevProc = (rank - 1 + size) % size;
  const int nextProc = (rank + 1) % size;

  std::string prevName = participantName + std::to_string(prevProc);
  std::string thisName = participantName + std::to_string(rank);
  std::string nextName = participantName + std::to_string(nextProc);
  if ((rank % 2) == 0) {
    left.prepareEstablishment(prevName, thisName);
    left.acceptConnection(prevName, thisName, tag, 0);
    left.cleanupEstablishment(prevName, thisName);

    right.requestConnection(thisName, nextName, tag, 0, 1);
  } else {
    right.requestConnection(thisName, nextName, tag, 0, 1);

    left.prepareEstablishment(prevName, thisName);
    left.acceptConnection(prevName, thisName, tag, 0);
    left.cleanupEstablishment(prevName, thisName);
  }
}

} // namespace precice::com
