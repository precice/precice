#ifndef PRECICE_NO_SOCKETS

#ifndef PRECICE_COM_SOCKET_REQUEST_HPP_
#define PRECICE_COM_SOCKET_REQUEST_HPP_

#include "Request.hpp"

#include <condition_variable>
#include <mutex>

namespace precice {
namespace com {
class SocketRequest : public Request {
public:
  SocketRequest();

  void complete();

  bool test();

  void wait();

private:
  bool _complete;

  std::condition_variable _completeCondition;
  std::mutex _completeMutex;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_SOCKET_REQUEST_HPP_ */

#endif // not PRECICE_NO_SOCKETS
