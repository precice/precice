#include "SocketRequest.hpp"

namespace precice::com {

void SocketRequest::complete()
{
  {
    std::lock_guard<std::mutex> lock(_completeMutex);

    _complete = true;
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

  // Lock is acquired when the predicate is evaluated.
  _completeCondition.wait(lock, [this] { return _complete; });
}
} // namespace precice::com
