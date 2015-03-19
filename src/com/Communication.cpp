// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#include "Communication.hpp"

#include "Request.hpp"

namespace precice {
namespace com {
void
Communication::broadcast() {
}

void
Communication::broadcast(int* itemsToSend, int size) {
  std::vector<com::PtrRequest> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (int rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(itemsToSend, size, rank + _rankOffset);

    requests.push_back(request);
  }

  for (auto request : requests) {
    request->wait();
  }
}

void
Communication::broadcast(int* itemsToReceive, int size, int rankBroadcaster) {
  receive(itemsToReceive, size, rankBroadcaster + _rankOffset);
}

void
Communication::broadcast(int itemToSend) {
  std::vector<com::PtrRequest> requests;

  requests.reserve(getRemoteCommunicatorSize());

  for (int rank = 0; rank < getRemoteCommunicatorSize(); ++rank) {
    auto request = aSend(&itemToSend, rank + _rankOffset);

    requests.push_back(request);
  }

  for (auto request : requests) {
    request->wait();
  }
}

void
Communication::broadcast(int& itemToReceive, int rankBroadcaster) {
  receive(itemToReceive, rankBroadcaster + _rankOffset);
}
}
} // namespace precice, com
