#include "SocketCommunication.hpp"
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <sstream>
#include "ConnectionInfoPublisher.hpp"
#include "SocketRequest.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace com {

namespace asio = boost::asio;

SocketCommunication::SocketCommunication(unsigned short     portNumber,
                                         bool               reuseAddress,
                                         std::string const &networkName,
                                         std::string const &addressDirectory)
    : _portNumber(portNumber),
      _reuseAddress(reuseAddress),
      _networkName(networkName),
      _addressDirectory(addressDirectory),
      _ioService(new IOService)
{
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

SocketCommunication::SocketCommunication(std::string const &addressDirectory)
    : SocketCommunication(0, false, "lo", addressDirectory)
{
}

SocketCommunication::~SocketCommunication()
{
  PRECICE_TRACE(_isConnected);
  closeConnection();
}

size_t SocketCommunication::getRemoteCommunicatorSize()
{
  PRECICE_TRACE();
  PRECICE_ASSERT(isConnected());
  return _sockets.size();
}

void SocketCommunication::acceptConnection(std::string const &acceptorName,
                                           std::string const &requesterName,
                                           std::string const &tag,
                                           int                acceptorRank,
                                           int                rankOffset)
{
  PRECICE_TRACE(acceptorName, requesterName);

  PRECICE_ASSERT(not isConnected());

  setRankOffset(rankOffset);

  std::string address;

  try {
    std::string ipAddress = getIpAddress();
    PRECICE_CHECK(not ipAddress.empty(), "Network \"" << _networkName << "\" not found for socket connection!");

    using asio::ip::tcp;

    tcp::acceptor acceptor(*_ioService);
    tcp::endpoint endpoint(tcp::v4(), _portNumber);

    acceptor.open(endpoint.protocol());
    acceptor.set_option(tcp::acceptor::reuse_address(_reuseAddress));
    acceptor.bind(endpoint);
    acceptor.listen();

    _portNumber = acceptor.local_endpoint().port();
    address     = ipAddress + ":" + std::to_string(_portNumber);
    ConnectionInfoWriter conInfo(acceptorName, requesterName, tag, _addressDirectory);
    conInfo.write(address);
    PRECICE_DEBUG("Accept connection at " << address);

    int peerCurrent               = 0;  // Current peer to connect to
    int peerCount                 = -1; // The total count of peers (initialized in the first iteration)
    int requesterCommunicatorSize = -1;

    do {
      auto socket = std::make_shared<Socket>(*_ioService);

      acceptor.accept(*socket);
      PRECICE_DEBUG("Accepted connection at " << address);
      _isConnected = true;

      int requesterRank = -1;

      asio::read(*socket, asio::buffer(&requesterRank, sizeof(int)));

      PRECICE_CHECK(_sockets.count(requesterRank) == 0,
                    "Duplicate request to connect by same rank (" << requesterRank << ")!");

      _sockets[requesterRank] = socket;
      send(acceptorRank, requesterRank);
      receive(requesterCommunicatorSize, requesterRank);

      // Initialize the count of peers to connect to
      if (peerCurrent == 0) {
        peerCount = requesterCommunicatorSize;
      }

      PRECICE_CHECK(requesterCommunicatorSize == peerCount,
                    "Requester communicator sizes are inconsistent!");
      PRECICE_CHECK(requesterCommunicatorSize > 0,
                    "Requester communicator size has to be > 0!");
    } while (++peerCurrent < requesterCommunicatorSize);

    acceptor.close();
  } catch (std::exception &e) {
    PRECICE_ERROR("Accepting connection at " << address << " failed: " << e.what());
  }

  // NOTE:
  // Keep IO service running so that it fires asynchronous handlers from another thread.
  _work   = std::make_shared<asio::io_service::work>(*_ioService);
  _thread = std::thread([this] { _ioService->run(); });
}

void SocketCommunication::acceptConnectionAsServer(std::string const &acceptorName,
                                                   std::string const &requesterName,
                                                   std::string const &tag,
                                                   int                acceptorRank,
                                                   int                requesterCommunicatorSize)
{
  PRECICE_TRACE(acceptorName, requesterName, acceptorRank, requesterCommunicatorSize);
  PRECICE_CHECK(requesterCommunicatorSize > 0, "Requester communicator size has to be > 0!");
  PRECICE_ASSERT(not isConnected());

  std::string address;

  try {
    std::string ipAddress = getIpAddress();
    PRECICE_CHECK(not ipAddress.empty(), "Network \"" << _networkName << "\" not found for socket connection!");

    using asio::ip::tcp;

    tcp::acceptor acceptor(*_ioService);
    {
      tcp::endpoint endpoint(tcp::v4(), _portNumber);

      acceptor.open(endpoint.protocol());
      acceptor.set_option(tcp::acceptor::reuse_address(_reuseAddress));
      acceptor.bind(endpoint);
      acceptor.listen();

      _portNumber = acceptor.local_endpoint().port();
    }

    address = ipAddress + ":" + std::to_string(_portNumber);
    ConnectionInfoWriter conInfo(acceptorName, requesterName, tag, acceptorRank, _addressDirectory);
    conInfo.write(address);

    PRECICE_DEBUG("Accepting connection at " << address);

    for (int connection = 0; connection < requesterCommunicatorSize; ++connection) {
      auto socket = std::make_shared<Socket>(*_ioService);
      acceptor.accept(*socket);
      PRECICE_DEBUG("Accepted connection at " << address);
      _isConnected = true;

      int requesterRank;
      asio::read(*socket, asio::buffer(&requesterRank, sizeof(int)));
      _sockets[requesterRank] = socket;
    }

    acceptor.close();
  } catch (std::exception &e) {
    PRECICE_ERROR("Accepting connection at " << address << " failed: " << e.what());
  }

  // NOTE: Keep IO service running so that it fires asynchronous handlers from another thread.
  _work   = std::make_shared<asio::io_service::work>(*_ioService);
  _thread = std::thread([this] { _ioService->run(); });
}

void SocketCommunication::requestConnection(std::string const &acceptorName,
                                            std::string const &requesterName,
                                            std::string const &tag,
                                            int                requesterRank,
                                            int                requesterCommunicatorSize)
{
  PRECICE_TRACE(acceptorName, requesterName);
  PRECICE_ASSERT(not isConnected());

  ConnectionInfoReader conInfo(acceptorName, requesterName, tag, _addressDirectory);
  std::string const    address = conInfo.read();
  PRECICE_DEBUG("Request connection to " << address);
  auto const        sepidx     = address.find(':');
  std::string const ipAddress  = address.substr(0, sepidx);
  std::string const portNumber = address.substr(sepidx + 1);
  _portNumber                  = static_cast<unsigned short>(std::stoul(portNumber));

  try {
    auto socket = std::make_shared<Socket>(*_ioService);

    using asio::ip::tcp;

    tcp::resolver::query query(tcp::v4(), ipAddress, portNumber, tcp::resolver::query::canonical_name);

    while (not isConnected()) {
      tcp::resolver                resolver(*_ioService);
      tcp::resolver::endpoint_type endpoint = *(resolver.resolve(query));
      boost::system::error_code    error    = asio::error::host_not_found;
      socket->connect(endpoint, error);

      _isConnected = not error;

      if (not isConnected()) {
        // Wait a little, since after a couple of ten-thousand trials the system
        // seems to get confused and the requester connects wrongly to itself.
        boost::asio::deadline_timer timer(*_ioService, boost::posix_time::milliseconds(1));
        timer.wait();
      }
    }

    PRECICE_DEBUG("Requested connection to " << address);

    asio::write(*socket, asio::buffer(&requesterRank, sizeof(int)));

    int acceptorRank = -1;
    asio::read(*socket, asio::buffer(&acceptorRank, sizeof(int)));
    _sockets[0] = socket; // should be acceptorRank instead of 0, likewise all communication below

    send(requesterCommunicatorSize, 0);

  } catch (std::exception &e) {
    PRECICE_ERROR("Requesting connection to " << address << " failed: " << e.what());
  }

  // NOTE: Keep IO service running so that it fires asynchronous handlers from another thread.
  _work   = std::make_shared<asio::io_service::work>(*_ioService);
  _thread = std::thread([this] { _ioService->run(); });
}

void SocketCommunication::requestConnectionAsClient(std::string const &  acceptorName,
                                                    std::string const &  requesterName,
                                                    std::string const &  tag,
                                                    std::set<int> const &acceptorRanks,
                                                    int                  requesterRank)

{
  PRECICE_TRACE(acceptorName, requesterName, acceptorRanks, requesterRank);
  PRECICE_ASSERT(not isConnected());

  for (auto const &acceptorRank : acceptorRanks) {
    _isConnected = false;
    ConnectionInfoReader conInfo(acceptorName, requesterName, tag, acceptorRank, _addressDirectory);
    std::string const    address    = conInfo.read();
    auto const           sepidx     = address.find(':');
    std::string const    ipAddress  = address.substr(0, sepidx);
    std::string const    portNumber = address.substr(sepidx + 1);
    _portNumber                     = static_cast<unsigned short>(std::stoul(portNumber));

    try {
      auto socket = std::make_shared<Socket>(*_ioService);

      using asio::ip::tcp;

      PRECICE_DEBUG("Requesting connection to " << ipAddress << ", port " << portNumber);

      tcp::resolver::query query(tcp::v4(), ipAddress, portNumber);

      while (not isConnected()) {
        tcp::resolver             resolver(*_ioService);
        tcp::resolver::iterator   endpoint_iterator = resolver.resolve(query);
        boost::system::error_code error             = asio::error::host_not_found;
        boost::asio::connect(*socket, endpoint_iterator, error);

        _isConnected = not error;

        if (not isConnected()) {
          // Wait a little, since after a couple of ten-thousand trials the system
          // seems to get confused and the requester connects wrongly to itself.
          boost::asio::deadline_timer timer(*_ioService, boost::posix_time::milliseconds(1));
          timer.wait();
        }
      }

      PRECICE_DEBUG("Requested connection to " << address << ", rank = " << acceptorRank);
      _sockets[acceptorRank] = socket;
      send(requesterRank, acceptorRank); // send my rank

    } catch (std::exception &e) {
      PRECICE_ERROR("Requesting connection to " << address << " failed: " << e.what());
    }
  }
  // NOTE: Keep IO service running so that it fires asynchronous handlers from another thread.
  _work   = std::make_shared<asio::io_service::work>(*_ioService);
  _thread = std::thread([this] { _ioService->run(); });
}

void SocketCommunication::closeConnection()
{
  PRECICE_TRACE();

  if (not isConnected())
    return;

  if (_thread.joinable()) {
    _work.reset();
    _ioService->stop();
    _thread.join();
  }

  for (auto &socket : _sockets) {
    PRECICE_ASSERT(socket.second->is_open());
    socket.second->shutdown(Socket::shutdown_both);
    socket.second->close();
  }

  _isConnected = false;
}

void SocketCommunication::send(std::string const &itemToSend, int rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  size_t size = itemToSend.size() + 1;
  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&size, sizeof(size_t)));
    asio::write(*_sockets[rankReceiver], asio::buffer(itemToSend.c_str(), size));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send failed: " << e.what());
  }
}

void SocketCommunication::send(const int *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size, rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(itemsToSend, size * sizeof(int)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send failed: " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const int *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size, rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  PtrRequest request(new SocketRequest);

  _queue.dispatch(_sockets[rankReceiver],
                  asio::buffer(itemsToSend, size * sizeof(int)),
                  [request] {
                    std::static_pointer_cast<SocketRequest>(request)->complete();
                  });
  return request;
}

void SocketCommunication::send(const double *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size, rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(itemsToSend, size * sizeof(double)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send failed: " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const double *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size, rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  PtrRequest request(new SocketRequest);

  _queue.dispatch(_sockets[rankReceiver],
                  asio::buffer(itemsToSend, size * sizeof(double)),
                  [request] {
                    std::static_pointer_cast<SocketRequest>(request)->complete();
                  });
  return request;
}

PtrRequest SocketCommunication::aSend(std::vector<double> const &itemsToSend, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  PtrRequest request(new SocketRequest);

  _queue.dispatch(_sockets[rankReceiver],
                  asio::buffer(itemsToSend),
                  [request] {
                    std::static_pointer_cast<SocketRequest>(request)->complete();
                  });
  return request;
}

void SocketCommunication::send(double itemToSend, int rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&itemToSend, sizeof(double)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send failed: " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const double &itemToSend, int rankReceiver)
{
  return aSend(&itemToSend, 1, rankReceiver);
}

void SocketCommunication::send(int itemToSend, int rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver)
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&itemToSend, sizeof(int)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send failed: " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const int &itemToSend, int rankReceiver)
{
  return aSend(&itemToSend, 1, rankReceiver);
}

void SocketCommunication::send(bool itemToSend, int rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&itemToSend, sizeof(bool)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send failed: " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const bool &itemToSend, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  PtrRequest request(new SocketRequest);

  _queue.dispatch(_sockets[rankReceiver],
                  asio::buffer(&itemToSend, sizeof(bool)),
                  [request] {
                    std::static_pointer_cast<SocketRequest>(request)->complete();
                  });
  return request;
}

void SocketCommunication::receive(std::string &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  size_t size = 0;

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&size, sizeof(size_t)));
    std::vector<char> msg(size);
    asio::read(*_sockets[rankSender], asio::buffer(msg.data(), size));
    itemToReceive = msg.data();
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }
}

void SocketCommunication::receive(int *itemsToReceive, int size, int rankSender)
{
  PRECICE_TRACE(size, rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(itemsToReceive, size * sizeof(int)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }
}

void SocketCommunication::receive(double *itemsToReceive, int size, int rankSender)
{
  PRECICE_TRACE(size, rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(itemsToReceive, size * sizeof(double)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }
}

PtrRequest SocketCommunication::aReceive(double *itemsToReceive,
                                         int     size,
                                         int     rankSender)
{
  PRECICE_TRACE(size, rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  PtrRequest request(new SocketRequest);

  try {
    asio::async_read(*_sockets[rankSender],
                     asio::buffer(itemsToReceive, size * sizeof(double)),
                     [request](boost::system::error_code const &, std::size_t) {
                       std::static_pointer_cast<SocketRequest>(request)->complete();
                     });
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }

  return request;
}

PtrRequest SocketCommunication::aReceive(std::vector<double> &itemsToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  PtrRequest request(new SocketRequest);

  try {
    asio::async_read(*_sockets[rankSender],
                     asio::buffer(itemsToReceive),
                     [request](boost::system::error_code const &, std::size_t) {
                       std::static_pointer_cast<SocketRequest>(request)->complete();
                     });
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }

  return request;
}

void SocketCommunication::receive(double &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&itemToReceive, sizeof(double)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }
}

PtrRequest SocketCommunication::aReceive(double &itemToReceive, int rankSender)
{
  return aReceive(&itemToReceive, 1, rankSender);
}

void SocketCommunication::receive(int &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&itemToReceive, sizeof(int)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }
}

PtrRequest SocketCommunication::aReceive(int &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT((rankSender >= 0) && (rankSender < (int) _sockets.size()),
                 rankSender, _sockets.size());
  PRECICE_ASSERT(isConnected());

  PtrRequest request(new SocketRequest);

  try {
    asio::async_read(*_sockets[rankSender],
                     asio::buffer(&itemToReceive, sizeof(int)),
                     [request](boost::system::error_code const &, std::size_t) {
                       std::static_pointer_cast<SocketRequest>(request)->complete();
                     });
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }

  return request;
}

void SocketCommunication::receive(bool &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&itemToReceive, sizeof(bool)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }
}

PtrRequest SocketCommunication::aReceive(bool &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  PtrRequest request(new SocketRequest);

  try {
    asio::async_read(*_sockets[rankSender],
                     asio::buffer(&itemToReceive, sizeof(bool)),
                     [request](boost::system::error_code const &, std::size_t) {
                       std::static_pointer_cast<SocketRequest>(request)->complete();
                     });
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }

  return request;
}

void SocketCommunication::send(std::vector<int> const &v, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  size_t size = v.size();
  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&size, sizeof(size_t)));
    asio::write(*_sockets[rankReceiver], asio::buffer(v));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send failed: " << e.what());
  }
}

void SocketCommunication::receive(std::vector<int> &v, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  size_t size = 0;

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&size, sizeof(size_t)));
    v.resize(size);
    asio::read(*_sockets[rankSender], asio::buffer(v));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }
}

void SocketCommunication::send(std::vector<double> const &v, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver);

  rankReceiver = rankReceiver - _rankOffset;

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  size_t size = v.size();
  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&size, sizeof(size_t)));
    asio::write(*_sockets[rankReceiver], asio::buffer(v));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send failed: " << e.what());
  }
}

void SocketCommunication::receive(std::vector<double> &v, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = rankSender - _rankOffset;

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  size_t size = 0;

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&size, sizeof(size_t)));
    v.resize(size);
    asio::read(*_sockets[rankSender], asio::buffer(v));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive failed: " << e.what());
  }
}

std::string SocketCommunication::getIpAddress()
{
  PRECICE_TRACE();
  std::ostringstream oss;

#ifdef _WIN32
  oss << "127.0.0.1";
#else
  int querySocket = socket(AF_INET, SOCK_STREAM, 0);

  PRECICE_ASSERT(querySocket >= 0);

  PRECICE_DEBUG("Looking for IP address of network \"" << _networkName << "\"");

  struct ifreq         request;
  struct if_nameindex *nameInterface = if_nameindex();

  PRECICE_ASSERT(nameInterface);

  struct if_nameindex *itNameInterface = nameInterface;

  while (itNameInterface && itNameInterface->if_name) {
    PRECICE_CHECK(strlen(itNameInterface->if_name) < IFNAMSIZ,
                  "Network interface \" " << itNameInterface->if_name << "\" has too long name");

    strncpy(request.ifr_name,
            itNameInterface->if_name,
            IFNAMSIZ); // Copy interface name

    int ifNameLength = strlen(itNameInterface->if_name);

    request.ifr_name[ifNameLength] = 0; // Add C-string 0

    if (ioctl(querySocket, SIOCGIFADDR, &request) >= 0) {
      PRECICE_DEBUG(itNameInterface->if_name
                    << ": " << inet_ntoa(((struct sockaddr_in *) &request.ifr_addr)->sin_addr));

      if (strcmp(itNameInterface->if_name, _networkName.c_str()) == 0) {
        oss << inet_ntoa(((struct sockaddr_in *) &request.ifr_addr)->sin_addr);
      }
    } else {
      PRECICE_CHECK(strcmp(itNameInterface->if_name, _networkName.c_str()) != 0,
                    "Could not obtain network IP from "
                        << "network \"" << itNameInterface->if_name << "\"");
    }

    itNameInterface++;
  }

  if_freenameindex(nameInterface);
  close(querySocket);
#endif

  return oss.str();
}

} // namespace com
} // namespace precice
