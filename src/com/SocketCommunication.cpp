#include "SocketCommunication.hpp"
#include <algorithm>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <sstream>
#include <stdexcept>
#include <utility>
#include "ConnectionInfoPublisher.hpp"
#include "SocketRequest.hpp"
#include "logging/LogMacros.hpp"
#include "utils/assertion.hpp"
#include "utils/networking.hpp"

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
    : SocketCommunication(0, false, utils::networking::loopbackInterfaceName(), addressDirectory)
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

      PRECICE_ASSERT(_sockets.count(requesterRank) == 0,
                     "Rank " << requesterRank << " has already been connected. Duplicate requests are not allowed.");

      _sockets[requesterRank] = socket;
      // send and receive expect a rank from the acceptor perspective.
      // Thus we need to apply given rankOffset before passing it to send/receive.
      // This is essentially the inverse of adjustRank().
      auto adjustedRequesterRank = requesterRank + rankOffset;
      send(acceptorRank, adjustedRequesterRank);
      receive(requesterCommunicatorSize, adjustedRequesterRank);

      // Initialize the count of peers to connect to
      if (peerCurrent == 0) {
        peerCount = requesterCommunicatorSize;
      }

      PRECICE_ASSERT(requesterCommunicatorSize > 0,
                     "Requester communicator size is " << requesterCommunicatorSize << " which is invalid.");
      PRECICE_ASSERT(requesterCommunicatorSize == peerCount,
                     "Current requester size from rank " << requesterRank << " is " << requesterCommunicatorSize << " but should be " << peerCount);
    } while (++peerCurrent < requesterCommunicatorSize);

    acceptor.close();
  } catch (std::exception &e) {
    PRECICE_ERROR("Accepting a socket connection at " << address << " failed with the system error: " << e.what());
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
  PRECICE_ASSERT(requesterCommunicatorSize >= 0, "Requester communicator size has to be positve.");
  PRECICE_ASSERT(not isConnected());

  if (requesterCommunicatorSize == 0) {
    PRECICE_DEBUG("Accepting no connections.");
    _isConnected = true;
    return;
  }

  std::string address;

  try {
    std::string ipAddress = getIpAddress();
    PRECICE_ASSERT(not ipAddress.empty(), "Network \"" << _networkName << "\" not found for socket connection!");

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
    PRECICE_ERROR("Accepting a socket connection at " << address << " failed with the system error: " << e.what());
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
    PRECICE_ERROR("Requesting a socket connection at " << address << " failed with the system error: " << e.what());
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
      PRECICE_ERROR("Requesting a socket connection at " << address << " failed with the system error: " << e.what());
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

    try {
      socket.second->shutdown(Socket::shutdown_send);
    } catch (std::exception &e) {
      PRECICE_WARN("Socket shutdown failed with system error: " << e.what());
    }
    socket.second->close();
  }

  _isConnected = false;
}

void SocketCommunication::send(std::string const &itemToSend, int rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);

  rankReceiver = adjustRank(rankReceiver);

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  size_t size = itemToSend.size() + 1;
  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&size, sizeof(size_t)));
    asio::write(*_sockets[rankReceiver], asio::buffer(itemToSend.c_str(), size));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send using sockets failed with system error: " << e.what());
  }
}

void SocketCommunication::send(const int *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size, rankReceiver);

  rankReceiver = adjustRank(rankReceiver);

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(itemsToSend, size * sizeof(int)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send using sockets failed with system error: " << e.what());
  }
}

void SocketCommunication::prepareEstablishment(std::string const &acceptorName,
                                               std::string const &requesterName)
{
  using namespace boost::filesystem;
  path dir = com::impl::localDirectory(acceptorName, requesterName, _addressDirectory);
  PRECICE_DEBUG("Creating connection exchange directory " << dir);
  try {
    create_directories(dir);
  } catch (const boost::filesystem::filesystem_error &e) {
    PRECICE_WARN("Creating directory for connection info failed with filesystem error: " << e.what());
  }
}

void SocketCommunication::cleanupEstablishment(std::string const &acceptorName,
                                               std::string const &requesterName)
{
  using namespace boost::filesystem;
  path dir = com::impl::localDirectory(acceptorName, requesterName, _addressDirectory);
  PRECICE_DEBUG("Removing connection exchange directory " << dir);
  try {
    remove_all(dir);
  } catch (const boost::filesystem::filesystem_error &e) {
    PRECICE_WARN("Cleaning up connection info failed with filesystem error " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const int *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size, rankReceiver);

  rankReceiver = adjustRank(rankReceiver);

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

  rankReceiver = adjustRank(rankReceiver);

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(itemsToSend, size * sizeof(double)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send using sockets failed with system error: " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const double *itemsToSend, int size, int rankReceiver)
{
  PRECICE_TRACE(size, rankReceiver);

  rankReceiver = adjustRank(rankReceiver);

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

  rankReceiver = adjustRank(rankReceiver);

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

  rankReceiver = adjustRank(rankReceiver);

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&itemToSend, sizeof(double)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send using sockets failed with system error: " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const double &itemToSend, int rankReceiver)
{
  return aSend(&itemToSend, 1, rankReceiver);
}

void SocketCommunication::send(int itemToSend, int rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);

  rankReceiver = adjustRank(rankReceiver);

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver)
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&itemToSend, sizeof(int)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send using sockets failed with system error: " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const int &itemToSend, int rankReceiver)
{
  return aSend(&itemToSend, 1, rankReceiver);
}

void SocketCommunication::send(bool itemToSend, int rankReceiver)
{
  PRECICE_TRACE(itemToSend, rankReceiver);

  rankReceiver = adjustRank(rankReceiver);

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&itemToSend, sizeof(bool)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send using sockets failed with system error: " << e.what());
  }
}

PtrRequest SocketCommunication::aSend(const bool &itemToSend, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver);

  rankReceiver = adjustRank(rankReceiver);

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

  rankSender = adjustRank(rankSender);

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  size_t size = 0;

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&size, sizeof(size_t)));
    std::vector<char> msg(size);
    asio::read(*_sockets[rankSender], asio::buffer(msg.data(), size));
    itemToReceive = msg.data();
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }
}

void SocketCommunication::receive(int *itemsToReceive, int size, int rankSender)
{
  PRECICE_TRACE(size, rankSender);

  rankSender = adjustRank(rankSender);

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(itemsToReceive, size * sizeof(int)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }
}

void SocketCommunication::receive(double *itemsToReceive, int size, int rankSender)
{
  PRECICE_TRACE(size, rankSender);

  rankSender = adjustRank(rankSender);

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(itemsToReceive, size * sizeof(double)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }
}

PtrRequest SocketCommunication::aReceive(double *itemsToReceive,
                                         int     size,
                                         int     rankSender)
{
  PRECICE_TRACE(size, rankSender);

  rankSender = adjustRank(rankSender);

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
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }

  return request;
}

PtrRequest SocketCommunication::aReceive(std::vector<double> &itemsToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = adjustRank(rankSender);

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
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }

  return request;
}

void SocketCommunication::receive(double &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = adjustRank(rankSender);

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&itemToReceive, sizeof(double)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }
}

PtrRequest SocketCommunication::aReceive(double &itemToReceive, int rankSender)
{
  return aReceive(&itemToReceive, 1, rankSender);
}

void SocketCommunication::receive(int &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = adjustRank(rankSender);

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&itemToReceive, sizeof(int)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }
}

PtrRequest SocketCommunication::aReceive(int &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = adjustRank(rankSender);

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
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }

  return request;
}

void SocketCommunication::receive(bool &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = adjustRank(rankSender);

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&itemToReceive, sizeof(bool)));
  } catch (std::exception &e) {
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }
}

PtrRequest SocketCommunication::aReceive(bool &itemToReceive, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = adjustRank(rankSender);

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
    PRECICE_ERROR("Receive using sockets failed with system error: " << e.what());
  }

  return request;
}

void SocketCommunication::send(std::vector<int> const &v, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver);

  rankReceiver = adjustRank(rankReceiver);

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  size_t size = v.size();
  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&size, sizeof(size_t)));
    asio::write(*_sockets[rankReceiver], asio::buffer(v));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send using sockets failed with system error: " << e.what());
  }
}

void SocketCommunication::receive(std::vector<int> &v, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = adjustRank(rankSender);

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  size_t size = 0;

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&size, sizeof(size_t)));
    v.resize(size);
    asio::read(*_sockets[rankSender], asio::buffer(v));
  } catch (std::exception &e) {
    PRECICE_ERROR("Recieve using sockets failed with system error: " << e.what());
  }
}

void SocketCommunication::send(std::vector<double> const &v, int rankReceiver)
{
  PRECICE_TRACE(rankReceiver);

  rankReceiver = adjustRank(rankReceiver);

  PRECICE_ASSERT(rankReceiver >= 0, rankReceiver);
  PRECICE_ASSERT(isConnected());

  size_t size = v.size();
  try {
    asio::write(*_sockets[rankReceiver], asio::buffer(&size, sizeof(size_t)));
    asio::write(*_sockets[rankReceiver], asio::buffer(v));
  } catch (std::exception &e) {
    PRECICE_ERROR("Send using sockets failed with system error: " << e.what());
  }
}

void SocketCommunication::receive(std::vector<double> &v, int rankSender)
{
  PRECICE_TRACE(rankSender);

  rankSender = adjustRank(rankSender);

  PRECICE_ASSERT(rankSender >= 0, rankSender);
  PRECICE_ASSERT(isConnected());

  size_t size = 0;

  try {
    asio::read(*_sockets[rankSender], asio::buffer(&size, sizeof(size_t)));
    v.resize(size);
    asio::read(*_sockets[rankSender], asio::buffer(v));
  } catch (std::exception &e) {
    PRECICE_ERROR("Recieve using sockets failed with system error: " << e.what());
  }
}

#ifndef _WIN32
namespace {
struct Interface {
  unsigned int index;
  std::string  name;
  std::string  address;
};

std::vector<Interface> detectInterfaces()
{
  std::vector<Interface> interfaces;

  // Collect interface indices and names
  struct if_nameindex *nameInterface = if_nameindex();
  for (struct if_nameindex *itNameInterface = nameInterface; itNameInterface->if_index != 0; ++itNameInterface) {
    Interface interface;
    interface.index = itNameInterface->if_index;
    interface.name  = itNameInterface->if_name;
    interfaces.emplace_back(std::move(interface));
  }
  if_freenameindex(nameInterface);

  // Resolve addresses for interfaces
  for (auto &interface : interfaces) {
    struct ifreq request;
    strncpy(request.ifr_name,
            interface.name.c_str(),
            IFNAMSIZ - 1); // Copy interface name

    auto socketfd = socket(AF_INET, SOCK_STREAM, 0);
    auto err      = ioctl(socketfd, SIOCGIFADDR, &request);
    close(socketfd);
    if (err) {
      continue;
    }

    const char *addr = inet_ntoa(((struct sockaddr_in *) &request.ifr_addr)->sin_addr);
    if (!addr) {
      continue;
    }
    interface.address = addr;
  }

  return interfaces;
}
} // namespace
#endif

std::string SocketCommunication::getIpAddress()
{
  PRECICE_TRACE();

#ifdef _WIN32
  return "127.0.0.1";
#else

  PRECICE_DEBUG("Looking for IP address of network \"" << _networkName << "\"");

  auto interfaces = detectInterfaces();

  auto pos = std::find_if(interfaces.begin(), interfaces.end(),
                          [&](Interface const &interface) { return interface.name == _networkName; });
  if (pos == interfaces.end()) {
    PRECICE_DEBUG("There  NOTHIGN");
    std::ostringstream err;
    err << "Cannot find network interface \"" << _networkName << "\". Available interfaces are: ";
    for (const auto &interface : interfaces) {
      err << interface.name << ' ';
    }
    err << " Please check \"network\" attribues in your configuration file.";
    PRECICE_ERROR(err.str());
  }

  // Unconnected interfaces don't have an IP address.
  PRECICE_CHECK(not pos->address.empty(), "The interface \"" << _networkName << "\" does not have an IP address. Please select another interface.");

  PRECICE_DEBUG("Detected network IP address of interface \"" << _networkName << "\":  " << pos->address << '.');
  return pos->address;
#endif
}

} // namespace com
} // namespace precice
