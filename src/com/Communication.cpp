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

  for (int i = 0; i < size; i++) {
    itemsToReceive[i] = itemsToSend[i];
  }

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

/**
 * @attention This method modifies the input buffer.
 */
void Communication::reduceSum(int &itemsToSend, int &itemsToReceive)
{
  TRACE();

  itemsToReceive = itemsToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(&itemsToSend, 1, rank + _rankOffset);
    request->wait();
    itemsToReceive += itemsToSend;
  }
}

void Communication::reduceSum(int &itemsToSend, int &itemsToReceive, int rankMaster)
{
  TRACE();

  auto request = aSend(&itemsToSend, 1, rankMaster);
  request->wait();
}

void Communication::allreduceSum()
{
  ERROR("Not implemented!");
  TRACE();
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(double *itemsToSend, double *itemsToReceive, int size)
{
  TRACE(size);

  for (int i = 0; i < size; i++) {
    itemsToReceive[i] = itemsToSend[i];
  }

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemsToSend, size, rank + _rankOffset);
    request->wait();
    for (int i = 0; i < size; i++) {
      itemsToReceive[i] += itemsToSend[i];
    }
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests;
  requests.reserve(getRemoteCommunicatorSize());
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

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(double &itemsToSend, double &itemsToReceive)
{
  TRACE();

  itemsToReceive = itemsToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(&itemsToSend, 1, rank + _rankOffset);
    request->wait();
    itemsToReceive += itemsToSend;
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests;
  requests.reserve(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemsToReceive, 1, rank + _rankOffset);
    requests.push_back(request);
  }
  Request::wait(requests);
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(double &itemsToSend, double &itemsToReceive, int rankMaster)
{
  TRACE();

  auto request = aSend(&itemsToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemsToReceive, 1, rankMaster + _rankOffset);
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(int &itemsToSend, int &itemsToReceive)
{
  TRACE();

  itemsToReceive = itemsToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(&itemsToSend, 1, rank + _rankOffset);
    request->wait();
    itemsToReceive += itemsToSend;
  }

  // send reduced result to all slaves
  std::vector<PtrRequest> requests;
  requests.reserve(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemsToReceive, 1, rank + _rankOffset);
    requests.push_back(request);
  }
  Request::wait(requests);
}

/**
 * @attention This method modifies the input buffer.
 */
void Communication::allreduceSum(int &itemsToSend, int &itemsToReceive, int rankMaster)
{
  TRACE();

  auto request = aSend(&itemsToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemsToReceive, 1, rankMaster + _rankOffset);
}

void Communication::broadcast()
{
  TRACE();
}

void Communication::broadcast(int *itemsToSend, int size)
{
  TRACE(size);

  std::vector<PtrRequest> requests;

  requests.reserve(getRemoteCommunicatorSize());

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

  std::vector<PtrRequest> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemToSend, rank + _rankOffset);

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

  std::vector<PtrRequest> requests;

  requests.reserve(getRemoteCommunicatorSize());

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

  std::vector<PtrRequest> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemToSend, rank + _rankOffset);

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
} // namespace com
} // namespace precice
