#include <algorithm>
#include <memory>
#include <ostream>

#include "Communication.hpp"
#include "Request.hpp"
#include "logging/LogMacros.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace com {

void Communication::connectMasterSlaves(std::string const &participantName,
                                        std::string const &tag,
                                        int                rank,
                                        int                size)
{
  if (size == 1)
    return;

  std::string masterName = participantName + "Master";
  std::string slaveName  = participantName + "Slave";

  constexpr Rank rankOffset = 1;
  int            slavesSize = size - rankOffset;
  if (rank == 0) {
    PRECICE_INFO("Connecting Master to {} Slaves", slavesSize);
    prepareEstablishment(masterName, slaveName);
    acceptConnection(masterName, slaveName, tag, rank, rankOffset);
    cleanupEstablishment(masterName, slaveName);
  } else {
    int slaveRank = rank - rankOffset;
    PRECICE_INFO("Connecting Slave #{} to Master", slaveRank);
    requestConnection(masterName, slaveName, tag, slaveRank, slavesSize);
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
  // receive local results from slaves
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request = aReceive(received, rank + _rankOffset);
    request->wait();
    for (size_t i = 0; i < itemsToReceive.size(); i++) {
      itemsToReceive[i] += received[i];
    }
  }
}

void Communication::reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank rankMaster)
{
  PRECICE_TRACE(itemsToSend.size(), itemsToReceive.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());

  auto request = aSend(itemsToSend, rankMaster);
  request->wait();
}

void Communication::reduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }
}

void Communication::reduceSum(int itemToSend, int &itemToReceive, Rank rankMaster)
{
  PRECICE_TRACE();

  auto request = aSend(itemToSend, rankMaster);
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

  // send reduced result to all slaves
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
void Communication::allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank rankMaster)
{
  PRECICE_TRACE(itemsToSend.size(), itemsToReceive.size());
  PRECICE_ASSERT(itemsToSend.size() == itemsToReceive.size());

  reduceSum(itemsToSend, itemsToReceive, rankMaster);
  // receive reduced data from master
  receive(itemsToReceive, rankMaster + _rankOffset);
}

void Communication::allreduceSum(double itemToSend, double &itemToReceive)
{
  PRECICE_TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request   = aSend(itemToReceive, rank + _rankOffset);
    requests[rank] = request;
  }
  Request::wait(requests);
}

void Communication::allreduceSum(double itemToSend, double &itemsToReceive, Rank rankMaster)
{
  PRECICE_TRACE();

  auto request = aSend(itemToSend, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(itemsToReceive, rankMaster + _rankOffset);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (Rank rank : remoteCommunicatorRanks()) {
    auto request   = aSend(itemToReceive, rank + _rankOffset);
    requests[rank] = request;
  }
  Request::wait(requests);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive, Rank rankMaster)
{
  PRECICE_TRACE();

  auto request = aSend(itemToSend, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(itemToReceive, rankMaster + _rankOffset);
}

void Communication::broadcast(precice::span<const int> itemsToSend)
{
  PRECICE_TRACE(itemsToSend.size());

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (Rank rank : remoteCommunicatorRanks()) {
    auto request   = aSend(itemsToSend, rank + _rankOffset);
    requests[rank] = request;
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
    auto request   = aSend(itemToSend, rank + _rankOffset);
    requests[rank] = request;
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
    auto request   = aSend(itemsToSend, rank + _rankOffset);
    requests[rank] = request;
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
    auto request   = aSend(itemToSend, rank + _rankOffset);
    requests[rank] = request;
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

int Communication::adjustRank(Rank rank) const
{
  return rank - _rankOffset;
}

} // namespace com
} // namespace precice
