#include "Communication.hpp"
#include <algorithm>
#include <memory>
#include <ostream>
#include "Request.hpp"
#include "logging/LogMacros.hpp"

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

  constexpr int rankOffset = 1;
  int           slavesSize = size - rankOffset;
  if (rank == 0) {
    PRECICE_INFO("Connecting Master to " << slavesSize << " Slaves");
    prepareEstablishment(masterName, slaveName);
    acceptConnection(masterName, slaveName, tag, rank, rankOffset);
    cleanupEstablishment(masterName, slaveName);
  } else {
    int slaveRank = rank - rankOffset;
    PRECICE_INFO("Connecting Slave #" << slaveRank << " to Master");
    requestConnection(masterName, slaveName, tag, slaveRank, slavesSize);
  }
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::reduceSum(double const *itemsToSend, double *itemsToReceive, int size)
{
  PRECICE_TRACE(size);

  std::copy(itemsToSend, itemsToSend + size, itemsToReceive);

  std::vector<double> received(size);
  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(received.data(), size, rank + _rankOffset);
    request->wait();
    for (int i = 0; i < size; i++) {
      itemsToReceive[i] += received[i];
    }
  }
}

void Communication::reduceSum(double const *itemsToSend, double *itemsToReceive, int size, int rankMaster)
{
  PRECICE_TRACE(size);

  auto request = aSend(itemsToSend, size, rankMaster);
  request->wait();
}

void Communication::reduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }
}

void Communication::reduceSum(int itemToSend, int &itemToReceive, int rankMaster)
{
  PRECICE_TRACE();

  auto request = aSend(&itemToSend, 1, rankMaster);
  request->wait();
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(double const *itemsToSend, double *itemsToReceive, int size)
{
  PRECICE_TRACE(size);

  std::copy(itemsToSend, itemsToSend + size, itemsToReceive);

  std::vector<double> received(size);
  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(received.data(), size, rank + _rankOffset);
    request->wait();
    for (int i = 0; i < size; i++) {
      itemsToReceive[i] += received[i];
    }
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request   = aSend(itemsToReceive, size, rank + _rankOffset);
    requests[rank] = request;
  }
  Request::wait(requests);
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(double const *itemsToSend, double *itemsToReceive, int size, int rankMaster)
{
  PRECICE_TRACE(size);

  auto request = aSend(itemsToSend, size, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(itemsToReceive, size, rankMaster + _rankOffset);
}

void Communication::allreduceSum(double itemToSend, double &itemToReceive)
{
  PRECICE_TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(&itemToSend, 1, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request   = aSend(&itemToReceive, 1, rank + _rankOffset);
    requests[rank] = request;
  }
  Request::wait(requests);
}

void Communication::allreduceSum(double itemToSend, double &itemsToReceive, int rankMaster)
{
  PRECICE_TRACE();

  auto request = aSend(&itemToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemsToReceive, 1, rankMaster + _rankOffset);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive)
{
  PRECICE_TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemToSend, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request   = aSend(&itemToReceive, 1, rank + _rankOffset);
    requests[rank] = request;
  }
  Request::wait(requests);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive, int rankMaster)
{
  PRECICE_TRACE();

  auto request = aSend(&itemToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemToReceive, 1, rankMaster + _rankOffset);
}

void Communication::broadcast(const int *itemsToSend, int size)
{
  PRECICE_TRACE(size);

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request   = aSend(itemsToSend, size, rank + _rankOffset);
    requests[rank] = request;
  }

  Request::wait(requests);
}

void Communication::broadcast(int *itemsToReceive, int size, int rankBroadcaster)
{
  PRECICE_TRACE(size);

  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(int itemToSend)
{
  PRECICE_TRACE();

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request   = aSend(itemToSend, rank + _rankOffset);
    requests[rank] = request;
  }

  Request::wait(requests);
}

void Communication::broadcast(int &itemToReceive, int rankBroadcaster)
{
  PRECICE_TRACE();
  receive(itemToReceive, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(const double *itemsToSend, int size)
{
  PRECICE_TRACE(size);

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request   = aSend(itemsToSend, size, rank + _rankOffset);
    requests[rank] = request;
  }

  Request::wait(requests);
}

void Communication::broadcast(double *itemsToReceive,
                              int     size,
                              int     rankBroadcaster)
{
  PRECICE_TRACE(size);
  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(double itemToSend)
{
  PRECICE_TRACE();

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request   = aSend(itemToSend, rank + _rankOffset);
    requests[rank] = request;
  }

  Request::wait(requests);
}

void Communication::broadcast(double &itemToReceive, int rankBroadcaster)
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

void Communication::broadcast(bool &itemToReceive, int rankBroadcaster)
{
  PRECICE_TRACE();
  int item;
  broadcast(item, rankBroadcaster);
  itemToReceive = item;
}

void Communication::broadcast(std::vector<int> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(const_cast<int *>(v.data()), v.size()); // make it send vector
}

void Communication::broadcast(std::vector<int> &v, int rankBroadcaster)
{
  v.clear();
  int size = 0;
  broadcast(size, rankBroadcaster);
  v.resize(size);
  broadcast(v.data(), size, rankBroadcaster);
}

void Communication::broadcast(std::vector<double> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(v.data(), v.size()); // make it send vector
}

void Communication::broadcast(std::vector<double> &v, int rankBroadcaster)
{
  v.clear();
  int size = 0;
  broadcast(size, rankBroadcaster);
  v.resize(size);
  broadcast(v.data(), size, rankBroadcaster);
}

int Communication::adjustRank(int rank) const
{
  return rank - _rankOffset;
}

} // namespace com
} // namespace precice
