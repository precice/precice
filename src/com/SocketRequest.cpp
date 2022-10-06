#include "SocketRequest.hpp"

namespace precice::com {
SocketRequest::SocketRequest()
    : _complete(false)
{
}

void SocketRequest::complete(boost::system::error_code ec)
{
  {
    std::lock_guard<std::mutex> lock(_completeMutex);

    _complete = true;
    _ec       = ec;
  }

  _completeCondition.notify_one();
}

bool SocketRequest::test()
{
  std::lock_guard<std::mutex> lock(_completeMutex);

  return _complete;
}

void SocketRequest::wait()
{
  std::unique_lock<std::mutex> lock(_completeMutex);

  _completeCondition.wait(lock, [this] { return _complete; });
}
} // namespace precice::com
