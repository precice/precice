#include <algorithm>
#include <boost/asio.hpp>
#include <iosfwd>
#include <new>
#include <utility>

#include "SocketSendQueue.hpp"
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace com {
namespace asio = boost::asio;

/// If items are left in the queue upon destruction, something went really wrong.
SocketSendQueue::~SocketSendQueue()
{
  PRECICE_ASSERT(_itemQueue.empty(), "The SocketSendQueue is not empty upon destruction. "
                                     "Make sure it always outlives all the requests pushed onto it.");
}

void SocketSendQueue::dispatch(std::shared_ptr<Socket>      sock,
                               boost::asio::const_buffers_1 data,
                               std::function<void()>        callback)
{
  _itemQueue.push_back({std::move(sock), std::move(data), callback});
  process(); // if queue was previously empty, start it now.
}

void SocketSendQueue::process()
{
  std::lock_guard<std::mutex> lock(_sendMutex);
  if (!_ready || _itemQueue.empty())
    return;
  auto item = _itemQueue.front();
  _itemQueue.pop_front();
  _ready = false;
  asio::async_write(*(item.sock),
                    item.data,
                    [item, this](boost::system::error_code const &, std::size_t) {
                      item.callback();
                      this->_ready = true;
                      this->process();
                    });
}

} // namespace com
} // namespace precice
