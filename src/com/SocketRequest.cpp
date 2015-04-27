// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_SOCKETS

#include "SocketRequest.hpp"

namespace precice {
namespace com {
SocketRequest::SocketRequest() : _complete(false) {
}

void
SocketRequest::complete() {
  {
    std::lock_guard<std::mutex> lock(_completeMutex);

    _complete = true;
  }

  _completeCondition.notify_one();
}

bool
SocketRequest::test() {
  std::lock_guard<std::mutex> lock(_completeMutex);

  return _complete;
}

void
SocketRequest::wait() {
  std::unique_lock<std::mutex> lock(_completeMutex);

  _completeCondition.wait(lock, [this] { return _complete; });
}
}
} // namespace precice, com

#endif // not PRECICE_NO_SOCKETS
