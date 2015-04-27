// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "Communication.hpp"

#include "Request.hpp"

#include "utils/Globals.hpp"

namespace precice {
namespace com {
tarch::logging::Log Communication::_log(
    "precice::com::Communication");



void
Communication::broadcast() {
  preciceTrace("broadcast()");
}

void
Communication::broadcast(int* itemsToSend, int size) {
  preciceTrace1("broadcast(int*)", size);

  std::vector<Request::SharedPointer> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (int rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToSend, size, rank + _rankOffset);

    requests.push_back(request);
  }

  Request::wait(requests);
}

void
Communication::broadcast(int* itemsToReceive, int size, int rankBroadcaster) {
  preciceTrace1("broadcast(int*)", size);

  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void
Communication::broadcast(int itemToSend) {
  preciceTrace("broadcast(int)");

  std::vector<Request::SharedPointer> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (int rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
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
  preciceTrace1("broadcast(double*)", size);

  std::vector<Request::SharedPointer> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (int rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToSend, size, rank + _rankOffset);

    requests.push_back(request);
  }

  Request::wait(requests);
}

void
Communication::broadcast(double* itemsToReceive,
                         int size,
                         int rankBroadcaster) {
  preciceTrace1("broadcast(double*)", size);

  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void
Communication::broadcast(double itemToSend) {
  preciceTrace("broadcast(double)");

  std::vector<Request::SharedPointer> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (int rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
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
