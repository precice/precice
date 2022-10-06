#pragma once

#include <condition_variable>
#include <mutex>
#include "Request.hpp"

namespace precice {
namespace com {
class SocketRequest : public Request {
public:
  SocketRequest();

  void complete(boost::system::error_code ec);

  bool test() override;

  void wait() override;

  boost::system::error_code errorCode() const noexcept override
  {
    return _ec;
  }

private:
  boost::system::error_code _ec;
  bool                      _complete;

  std::condition_variable _completeCondition;
  std::mutex              _completeMutex;
};
} // namespace com
} // namespace precice
