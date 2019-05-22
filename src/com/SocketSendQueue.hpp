#pragma once

#include <deque>
#include <boost/asio.hpp>
#include <functional>
#include <mutex>
#include "logging/Logger.hpp"

namespace precice
{
namespace com
{

/// This Queue is intended for SocketCommunication to push requests which should be sent onto it.
/// It ensures that the invocations of asio::aSend are done serially.
class SocketSendQueue
{
public:

  using TCP           = boost::asio::ip::tcp;
  using SocketService = boost::asio::stream_socket_service<TCP>;
  using Socket        = boost::asio::basic_stream_socket<TCP, SocketService>;

  SocketSendQueue() = default;
  ~SocketSendQueue();

  SocketSendQueue(SocketSendQueue const &) =delete;
  SocketSendQueue& operator=(SocketSendQueue const &) =delete;

  /// Put data in the queue, start processing the queue.
  void dispatch(std::shared_ptr<Socket> sock, boost::asio::const_buffers_1 data, std::function<void()> callback);

private:
  logging::Logger _log{"com::SocketSendQueue"};

  /// This method can be called arbitrarily many times, but enough times to ensure the queue makes progress.
  void process();

  struct SendItem
  {
    std::shared_ptr<Socket> sock;
    boost::asio::const_buffers_1 data;
    std::function<void()> callback;
  };

  std::deque<SendItem> _itemQueue;
  std::mutex _sendMutex;
  bool _ready = true;
};

}}
