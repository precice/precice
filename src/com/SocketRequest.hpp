#pragma once

#include <condition_variable>
#include <mutex>
#include "Request.hpp"

namespace precice {
namespace com {
class SocketRequest : public Request {
public:
  SocketRequest();

  void complete();

  bool test() override;

  void wait() override;

private:
  bool _complete;

  std::condition_variable _completeCondition;
  std::mutex              _completeMutex;
};
} // namespace com
} // namespace precice
