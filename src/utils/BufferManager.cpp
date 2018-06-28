#include "BufferManager.hpp"
#include "assertion.hpp"

namespace precice {
namespace utils {


class release_visitor : public boost::static_visitor<>
{
public:
  template<typename T>
  void operator()(T & operand) const
  {
    operand.release();
  }
};

class is_valid_ptr_visitor : public boost::static_visitor<bool>
{
public:
  template<typename T>
  bool operator()(T & operand) const
  {
    return static_cast<bool>(operand);
  }
};

BufferManager & BufferManager::instance()
{
  static BufferManager instance;
  return instance;
}

void BufferManager::run()
{
  stopFlag = false;
  if (not thread.joinable()) {
    INFO("Buffer manager thread starting");
    thread = std::thread(&BufferManager::check, this);
  }
}

void BufferManager::stop()
{
  stopFlag = true;
  thread.join();
}

void BufferManager::wait()
{
  while (thread.joinable()) {
    if (bufferedRequests.empty())
      break;
    std::this_thread::yield();    
  }
}

void BufferManager::check()
{
  while (not stopFlag) {
    std::lock_guard<std::mutex> lock(bufferedRequestsMutex); // maybe too coarse locking
    for (auto it = bufferedRequests.begin(); it != bufferedRequests.end();) {
      // if (boost::apply_visitor(is_valid_ptr_visitor(), it->second) and it->first->test()) {
      if (it->first->test()) {
        // boost::apply_visitor(release_visitor(), it->second);
        it = bufferedRequests.erase(it);
      }
      else {
        ++it;
      }
    }
    std::this_thread::yield();
  }
}

}
}
