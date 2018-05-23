#include "Communication.hpp"
#include "Request.hpp"

namespace precice
{
namespace com
{
/**
 * @attention This method modifies the input buffer.
 */
void Communication::reduceSum(double *itemsToSend, double *itemsToReceive, int size)
{
  TRACE(size);

  std::copy(itemsToSend, itemsToSend + size, itemsToReceive);
  
  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemsToSend, size, rank + _rankOffset);
    request->wait();
    for (int i = 0; i < size; i++) {
      itemsToReceive[i] += itemsToSend[i];
    }
  }
}

void Communication::reduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster)
{
  TRACE(size);

  auto request = aSend(itemsToSend, size, rankMaster);
  request->wait();
}

void Communication::reduceSum(int itemToSend, int &itemToReceive)
{
  TRACE();

  itemToReceive = itemToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(&itemToSend, 1, rank + _rankOffset);
    request->wait();
    itemToReceive += itemToSend;
  }
}

void Communication::reduceSum(int itemToSend, int &itemToReceive, int rankMaster)
{
  TRACE();

  auto request = aSend(&itemToSend, 1, rankMaster);
  request->wait();
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(double *itemsToSend, double *itemsToReceive, int size)
{
  TRACE(size);

  std::copy(itemsToSend, itemsToSend + size, itemsToReceive);
  
  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemsToSend, size, rank + _rankOffset);
    request->wait();
    for (int i = 0; i < size; i++) {
      itemsToReceive[i] += itemsToSend[i];
    }
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToReceive, size, rank + _rankOffset);
    requests.push_back(request);
  }
  Request::wait(requests);
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster)
{
  TRACE(size);

  auto request = aSend(itemsToSend, size, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(itemsToReceive, size, rankMaster + _rankOffset);
}

void Communication::allreduceSum(double itemToSend, double &itemToReceive)
{
  TRACE();

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
    auto request = aSend(&itemToReceive, 1, rank + _rankOffset);
    requests.push_back(request);
  }
  Request::wait(requests);
}

void Communication::allreduceSum(double itemToSend, double &itemsToReceive, int rankMaster)
{
  TRACE();

  auto request = aSend(&itemToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemsToReceive, 1, rankMaster + _rankOffset);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive)
{
  TRACE();

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
    auto request = aSend(&itemToReceive, 1, rank + _rankOffset);
    requests.push_back(request);
  }
  Request::wait(requests);
}

void Communication::allreduceSum(int itemToSend, int &itemToReceive, int rankMaster)
{
  TRACE();

  auto request = aSend(&itemToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemToReceive, 1, rankMaster + _rankOffset);
}

void Communication::broadcast(int *itemsToSend, int size)
{
  TRACE(size);

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToSend, size, rank + _rankOffset);

    requests.push_back(request);
  }

  Request::wait(requests);
}

void Communication::broadcast(int *itemsToReceive, int size, int rankBroadcaster)
{
  TRACE(size);

  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(int itemToSend)
{
  TRACE();

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemToSend, rank + _rankOffset);
    requests.push_back(request);
  }

  Request::wait(requests);
}

void Communication::broadcast(int &itemToReceive, int rankBroadcaster)
{
  TRACE();
  receive(itemToReceive, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(double *itemsToSend, int size)
{
  TRACE(size);

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToSend, size, rank + _rankOffset);

    requests.push_back(request);
  }

  Request::wait(requests);
}

void Communication::broadcast(double *itemsToReceive,
                              int     size,
                              int     rankBroadcaster)
{
  TRACE(size);
  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(double itemToSend)
{
  TRACE();

  std::vector<PtrRequest> requests(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemToSend, rank + _rankOffset);

    requests.push_back(request);
  }

  Request::wait(requests);
}

void Communication::broadcast(double &itemToReceive, int rankBroadcaster)
{
  TRACE();
  receive(itemToReceive, rankBroadcaster + _rankOffset);
}

void Communication::broadcast(bool itemToSend)
{
  TRACE();
  int item = itemToSend;
  broadcast(item);
}

void Communication::broadcast(bool &itemToReceive, int rankBroadcaster)
{
  TRACE();
  int item;
  broadcast(item, rankBroadcaster);
  itemToReceive = item;
}

void Communication::broadcast(std::vector<int> const &v)
{
  broadcast(static_cast<int>(v.size()));
  broadcast(const_cast<int*>(v.data()), v.size()); // make it send vector
}

void Communication::broadcast(std::vector<int> &v, int rankBroadcaster)
{
  v.clear();
  int size = 0;
  broadcast(size, rankBroadcaster);
  v.resize(size);
  broadcast(v.data(), size, rankBroadcaster);
}

/// @todo Reimplement more efficiently for MPI, e.g. using MPI_Probe
void Communication::send(std::vector<int> const &v, int rankReceiver)
{
  send(static_cast<int>(v.size()), rankReceiver);
  send(const_cast<int*>(v.data()), v.size(), rankReceiver);
}

/// @todo Reimplement more efficiently for MPI, e.g. using MPI_Probe
void Communication::receive(std::vector<int> &v, int rankSender)
{
  v.clear();
  int size = -1;
  receive(size, rankSender);
  v.resize(size);
  receive(v.data(), size, rankSender);
}


} // namespace com
} // namespace precice
