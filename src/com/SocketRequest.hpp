#pragma once
#ifndef PRECICE_NO_SOCKETS

#include "Request.hpp"

#include <condition_variable>
#include <mutex>

namespace precice
{
namespace com
{
class SocketRequest : public Request
{
public:
  SocketRequest();

  void complete();

  bool test();

  void wait();

private:
  bool _complete;

  std::condition_variable _completeCondition;
  std::mutex              _completeMutex;
};
}
} // namespace precice, com

#endif // not PRECICE_NO_SOCKETS
