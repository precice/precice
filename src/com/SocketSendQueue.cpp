#include "SocketSendQueue.hpp"

namespace precice
{
namespace com
{
namespace asio = boost::asio;

void SendQueue::push(std::shared_ptr<Socket> sock, boost::asio::const_buffers_1 data, std::function<void()> callback)
{
  _itemQueue.push_back({sock, data, callback});
  process(); //If queue was previously empty, start it now.
}
void SendQueue::process()
{
  std::lock_guard<std::mutex> lock(_queueMutex);
  if (!ready || _itemQueue.empty())
    return;
  auto item = _itemQueue.front();
  _itemQueue.pop_front();
  ready = false;
  asio::async_write(*(item.sock),
                    item.data,
                    [item, this](boost::system::error_code const &, std::size_t) {
                      item.callback();
                      this->ready = true;
                      this->process();
                    });
}
}}
