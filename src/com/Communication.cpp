// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "Communication.hpp"

#include "Request.hpp"

#include "utils/Globals.hpp"

namespace precice {
namespace com {
logging::Logger Communication::_log(
    "precice::com::Communication");



/**
 * Attention: this method modifies the input buffer.
 */
void
Communication::reduceSum(double* itemsToSend, double* itemsToReceive, int size) {
  preciceTrace("broadcast(double*)", size);

  for(int i = 0; i < size; i++){
    itemsToReceive[i] = itemsToSend[i];
  }

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemsToSend, size, rank + _rankOffset);
    request->wait();
    for(int i = 0; i < size; i++){
      itemsToReceive[i] += itemsToSend[i];
    }
  }
}


void
Communication::reduceSum(double* itemsToSend, double* itemsToReceive, int size, int rankMaster) {
  preciceTrace("allreduce(double*)", size);

  auto request = aSend(itemsToSend, size, rankMaster);
  request->wait();
}

/**
 * Attention: this method modifies the input buffer.
 */
void
Communication::reduceSum(int& itemsToSend, int& itemsToReceive) {
  preciceTrace("broadcast(int)");

  itemsToReceive = itemsToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(&itemsToSend, 1, rank + _rankOffset);
    request->wait();
    itemsToReceive += itemsToSend;
  }
}


void
Communication::reduceSum(int& itemsToSend, int& itemsToReceive, int rankMaster) {
  preciceTrace("allreduce(int)");

  auto request = aSend(&itemsToSend, 1, rankMaster);
  request->wait();
}

void
Communication::allreduceSum() {
  preciceTrace("allreduceSum()");
}

/**
 * Attention: this method modifies the input buffer.
 */
void
Communication::allreduceSum(double* itemsToSend, double* itemsToReceive, int size) {
  preciceTrace("broadcast(double*)", size);

  for(int i = 0; i < size; i++){
    itemsToReceive[i] = itemsToSend[i];
  }

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(itemsToSend, size, rank + _rankOffset);
    request->wait();
    for(int i = 0; i < size; i++){
      itemsToReceive[i] += itemsToSend[i];
    }
  }

  // send reduced result to all slaves
  std::vector<Request::SharedPointer> requests;
  requests.reserve(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToReceive, size, rank + _rankOffset);
    requests.push_back(request);
  }
  Request::wait(requests);
}

/**
 * Attention: this method modifies the input buffer.
 */
void
Communication::allreduceSum(double* itemsToSend, double* itemsToReceive, int size, int rankMaster) {
  preciceTrace("allreduce(double*)", size);

  auto request = aSend(itemsToSend, size, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(itemsToReceive, size, rankMaster + _rankOffset);
}

/**
 * Attention: this method modifies the input buffer.
 */
void
Communication::allreduceSum(double& itemsToSend, double& itemsToReceive) {
  preciceTrace("broadcast(double)");

  itemsToReceive = itemsToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(&itemsToSend, 1, rank + _rankOffset);
    request->wait();
    itemsToReceive += itemsToSend;
  }

  // send reduced result to all slaves
  std::vector<Request::SharedPointer> requests;
  requests.reserve(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemsToReceive, 1, rank + _rankOffset);
    requests.push_back(request);
  }
  Request::wait(requests);
}

/**
 * Attention: this method modifies the input buffer.
 */
void
Communication::allreduceSum(double& itemsToSend, double& itemsToReceive, int rankMaster) {
  preciceTrace("allreduce(double)");

  auto request = aSend(&itemsToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemsToReceive, 1, rankMaster + _rankOffset);
}

/**
 * Attention: this method modifies the input buffer.
 */
void
Communication::allreduceSum(int& itemsToSend, int& itemsToReceive) {
  preciceTrace("broadcast(int)");

  itemsToReceive = itemsToSend;

  // receive local results from slaves
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aReceive(&itemsToSend, 1, rank + _rankOffset);
    request->wait();
    itemsToReceive += itemsToSend;
  }

  // send reduced result to all slaves
  std::vector<Request::SharedPointer> requests;
  requests.reserve(getRemoteCommunicatorSize());
  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemsToReceive, 1, rank + _rankOffset);
    requests.push_back(request);
  }
  Request::wait(requests);
}

/**
 * Attention: this method modifies the input buffer.
 */
void
Communication::allreduceSum(int& itemsToSend, int& itemsToReceive, int rankMaster) {
  preciceTrace("allreduce(int)");

  auto request = aSend(&itemsToSend, 1, rankMaster);
  request->wait();
  // receive reduced data from master
  receive(&itemsToReceive, 1, rankMaster + _rankOffset);
}



void
Communication::broadcast() {
  preciceTrace("broadcast()");
}

void
Communication::broadcast(int* itemsToSend, int size) {
  preciceTrace("broadcast(int*)", size);

  std::vector<Request::SharedPointer> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToSend, size, rank + _rankOffset);

    requests.push_back(request);
  }

  Request::wait(requests);
}

void
Communication::broadcast(int* itemsToReceive, int size, int rankBroadcaster) {
  preciceTrace("broadcast(int*)", size);

  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void
Communication::broadcast(int itemToSend) {
  preciceTrace("broadcast(int)");

  std::vector<Request::SharedPointer> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemToSend, rank + _rankOffset);

    requests.push_back(request);
  }

  Request::wait(requests);
}

void
Communication::broadcast(int& itemToReceive, int rankBroadcaster) {
  preciceTrace("broadcast(int&)");

  receive(itemToReceive, rankBroadcaster + _rankOffset);
}

void
Communication::broadcast(double* itemsToSend, int size) {
  preciceTrace("broadcast(double*)", size);

  std::vector<Request::SharedPointer> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToSend, size, rank + _rankOffset);

    requests.push_back(request);
  }

  Request::wait(requests);
}

void
Communication::broadcast(double* itemsToReceive,
                         int size,
                         int rankBroadcaster) {
  preciceTrace("broadcast(double*)", size);
  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void
Communication::broadcast(double itemToSend) {
  preciceTrace("broadcast(double)");

  std::vector<Request::SharedPointer> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (size_t rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemToSend, rank + _rankOffset);

    requests.push_back(request);
  }

  Request::wait(requests);
}

void
Communication::broadcast(double& itemToReceive, int rankBroadcaster) {
  preciceTrace("broadcast(double&)");
  receive(itemToReceive, rankBroadcaster + _rankOffset);
}

void
Communication::broadcast(bool itemToSend) {
  preciceTrace("broadcast(bool)");

  int item = itemToSend;

  broadcast(item);
}

void
Communication::broadcast(bool& itemToReceive, int rankBroadcaster) {
  preciceTrace("broadcast(bool&)");

  int item;

  broadcast(item, rankBroadcaster);

  itemToReceive = item;
}
}
} // namespace precice, com
