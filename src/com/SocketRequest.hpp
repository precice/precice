#pragma once

#include <condition_variable>
#include <mutex>
#include "Request.hpp"

namespace precice::com {
class SocketRequest : public Request {
public:
  void complete();

  bool test() override;

  void wait() override;

private:
  bool _complete{false};

  std::condition_variable _completeCondition;
  std::mutex              _completeMutex;
};
} // namespace precice::com
