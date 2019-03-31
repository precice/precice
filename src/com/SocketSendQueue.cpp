#include "SocketSendQueue.hpp"

namespace precice
{
namespace com
{
namespace asio = boost::asio;

void SendQueue::dispatch(std::shared_ptr<Socket> sock,
        boost::asio::const_buffers_1 data,
        std::function<void()> callback)
{
  _itemQueue.push_back({std::move(sock), std::move(data), callback});
  process(); //If queue was previously empty, start it now.
}

/// This method can be called arbitrarily many times,
/// but enough times to ensure the queue makes progress.
void SendQueue::process()
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

/// If Items are left in the queue upon destruction, something went really wrong.
SendQueue::~SendQueue()
{
    if (not _itemQueue.empty())
        ERROR("The SocketSendQueue is not empty upon destruction.\n"
              "Make sure it always outlives all the requests pushed onto it");
}
}}
