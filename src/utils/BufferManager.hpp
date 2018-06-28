#pragma once

#include <memory>
#include <thread>
#include <list>
#include <mutex>
#include <boost/variant.hpp>
#include "com/Request.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace utils {


class BufferManager
{
public:
  /// Deleted copy operator for singleton pattern
  BufferManager(BufferManager const &) = delete;
  
  /// Deleted assigment operator for singleton pattern
  void operator=(BufferManager const &) = delete;

  static BufferManager & instance();

  /// Puts a new (Request, pointer) into the queue
  template<typename T>
  void put(std::shared_ptr<com::Request> request, std::shared_ptr<T> buffer)
  {
    TRACE();
    std::lock_guard<std::mutex> lock(bufferedRequestsMutex);
    bufferedRequests.emplace_back(request, buffer);
  }
  
  /// Starts the checking thread
  void run();

  /// Stops the checking thread
  void stop();

  /// Blocks until the queue is empty
  void wait();
  
  void check();

private:
  logging::Logger _log{"utils::BufferManager"};
  
  /// Private, empty constructor for singleton pattern
  BufferManager() {}

  using ptr_types = boost::variant< std::unique_ptr<std::vector<double>> >;

  // std::list<std::pair<std::unique_ptr<com::Request>, ptr_types>> bufferedRequests;
  std::list<std::pair<std::shared_ptr<com::Request>,
                      std::shared_ptr<std::vector<double>>>> bufferedRequests;
  
  std::thread thread;
  std::mutex bufferedRequestsMutex;
  bool stopFlag = false;
};


}
}
