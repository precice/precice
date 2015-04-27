// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_SOCKETS

#ifndef PRECICE_COM_SOCKET_REQUEST_HPP_
#define PRECICE_COM_SOCKET_REQUEST_HPP_

#include "Request.hpp"

#include <condition_variable>
#include <mutex>

namespace precice {
namespace com {
class SocketRequest : public Request {
public:
  SocketRequest();

  void complete();

  bool test();

  void wait();

private:
  bool _complete;

  std::condition_variable _completeCondition;
  std::mutex _completeMutex;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_SOCKET_REQUEST_HPP_ */

#endif // not PRECICE_NO_SOCKETS
