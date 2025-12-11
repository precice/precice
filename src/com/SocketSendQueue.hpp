#pragma once

#include <boost/asio.hpp>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include "logging/Logger.hpp"

namespace precice::com {

/// This Queue is intended for SocketCommunication to push requests which should be sent onto it.
/// It ensures that the invocations of asio::aSend are done serially.
class SocketSendQueue {
public:
  using Socket = boost::asio::ip::tcp::socket;

  SocketSendQueue() = default;
  ~SocketSendQueue();

  SocketSendQueue(SocketSendQueue const &)            = delete;
  SocketSendQueue &operator=(SocketSendQueue const &) = delete;

  /// Put data in the queue, start processing the queue.
  void dispatch(std::shared_ptr<Socket> sock, boost::asio::const_buffer data, std::function<void()> callback);

  /// Notifies the queue that the last asynchronous send operation has completed.
  void sendCompleted();

private:
  /// This method can be called arbitrarily many times, but enough times to ensure the queue makes progress.
  void process();

  struct SendItem {
    std::shared_ptr<Socket>   sock;
    boost::asio::const_buffer data;
    std::function<void()>     callback;
  };

  /// The queue, containing items to asynchronously send using boost.asio.
  std::deque<SendItem> _itemQueue;
  /// The mutex protecting access to the queue
  std::mutex _queueMutex{};
  /// Is the queue allowed to start another asynchronous send?
  bool _ready = true;
};

} // namespace precice::com
