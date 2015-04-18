// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NO_SOCKETS

#include "SocketCommunication.hpp"

#include "SocketRequest.hpp"

#include "utils/Publisher.hpp"

#include <boost/asio.hpp>
#include <boost/bind.hpp>

#include <sstream>

using precice::utils::Publisher;
using precice::utils::ScopedPublisher;

namespace precice {
namespace com {

namespace asio = boost::asio;

tarch::logging::Log SocketCommunication::_log(
    "precice::com::SocketCommunication");

SocketCommunication::SocketCommunication(unsigned short portNumber,
                                         bool reuseAddress,
                                         std::string const& networkName,
                                         std::string const& addressDirectory)
    : _portNumber(portNumber)
    , _reuseAddress(reuseAddress)
    , _networkName(networkName)
    , _addressDirectory(addressDirectory)
    , _isConnected(false)
    , _isClient(false)
    , _remoteCommunicatorSize(0)
    , _ioService(new IOService)
    , _sockets()
    , _queryWork()
    , _queryThread()
    , _clientQueries()
    , _clientQueryBuffers() {
  if (_addressDirectory.empty()) {
    _addressDirectory = ".";
  }
}

SocketCommunication::SocketCommunication(std::string const& addressDirectory)
    : SocketCommunication(0, false, "lo", addressDirectory) {
}

SocketCommunication::~SocketCommunication() {
  preciceTrace1("~SocketCommunication()", _isConnected);

  closeConnection();
}

bool
SocketCommunication::isConnected() {
  return _isConnected;
}

int
SocketCommunication::getRemoteCommunicatorSize() {
  preciceTrace("getRemoteCommunicatorSize()");

  assertion(isConnected());

  return _remoteCommunicatorSize;
}

void
SocketCommunication::acceptConnection(std::string const& nameAcceptor,
                                      std::string const& nameRequester,
                                      int acceptorProcessRank,
                                      int acceptorCommunicatorSize) {
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);
  preciceCheck(acceptorCommunicatorSize == 1,
               "acceptConnection()",
               "Acceptor of socket connection can only have one process!");

  _rank = acceptorProcessRank;

  using asio::ip::tcp;

  try {
    std::string ipAddress = getIpAddress();

    preciceCheck(not ipAddress.empty(),
                 "acceptConnection()",
                 "Network \"" << _networkName
                              << "\" not found for socket connection!");

    tcp::acceptor acceptor(*_ioService);

    {
      tcp::endpoint endpoint(tcp::v4(), _portNumber);

      acceptor.open(endpoint.protocol());
      acceptor.set_option(tcp::acceptor::reuse_address(_reuseAddress));
      acceptor.bind(endpoint);
      acceptor.listen();

      _portNumber = acceptor.local_endpoint().port();
    }

    std::string address(ipAddress + ":" + std::to_string(_portNumber));
    std::string addressFileName("." + nameRequester + "-" + nameAcceptor +
                                ".address");

    Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);

    ScopedPublisher p(addressFileName);

    p.write(address);

    // std::cout << address << std::endl;

    preciceDebug("Accept connection at " << address);

    PtrSocket socket(new tcp::socket(*_ioService));

    acceptor.accept(*socket);

    _isConnected = true;

    int remoteRank = -1;
    int remoteSize = 0;

    asio::read(*socket, asio::buffer((void*)&remoteRank, sizeof(int)));
    asio::read(*socket, asio::buffer((void*)&remoteSize, sizeof(int)));

    preciceCheck(remoteSize > 0,
                 "acceptConnection()",
                 "Requester communicator "
                     << "size has to be > 0!");

    _remoteCommunicatorSize = remoteSize;

    preciceDebug("Received rank=" << remoteRank << ", size=" << remoteSize);

    _sockets.resize(_remoteCommunicatorSize);
    _clientQueryBuffers.resize(_remoteCommunicatorSize);

    _sockets[remoteRank] = socket;

    send(acceptorProcessRank, remoteRank);
    send(acceptorCommunicatorSize, remoteRank);

    for (int i = 1; i < _remoteCommunicatorSize; ++i) {
      socket = PtrSocket(new Socket(*_ioService));
      acceptor.accept(*socket); // Waits until connection
      asio::read(*socket, asio::buffer((void*)&remoteRank, sizeof(int)));
      asio::read(*socket, asio::buffer((void*)&remoteSize, sizeof(int)));
      preciceCheck(remoteSize == _remoteCommunicatorSize,
                   "acceptConnection()",
                   "Remote communicator sizes are inconsistent!");
      preciceCheck(_sockets[remoteRank].use_count() == 0,
                   "acceptConnection()",
                   "Duplicate request to connect by same rank (" << remoteRank
                                                                 << ")!");
      _sockets[remoteRank] = socket;
      send(acceptorProcessRank, remoteRank);
      send(acceptorCommunicatorSize, remoteRank);
    }

    acceptor.close();
  } catch (std::exception& e) {
    preciceError("acceptConnection()",
                 "Accepting connection at port " << _portNumber
                                                 << " failed: " << e.what());
  }

  _queryWork = PtrWork(new asio::io_service::work(*_ioService));
  _queryThread = std::thread(&SocketCommunication::onThreadRun, this);
}

void
SocketCommunication::acceptConnectionAsServer(std::string const& nameAcceptor,
                                              std::string const& nameRequester,
                                              int requesterCommunicatorSize) {
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);

  _rank = 0;

  using asio::ip::tcp;

  try {
    std::string ipAddress = getIpAddress();

    preciceCheck(not ipAddress.empty(),
                 "acceptConnection()",
                 "Network \"" << _networkName
                              << "\" not found for socket connection!");

    tcp::acceptor acceptor(*_ioService);

    {
      tcp::endpoint endpoint(tcp::v4(), _portNumber);

      acceptor.open(endpoint.protocol());
      acceptor.set_option(tcp::acceptor::reuse_address(_reuseAddress));
      acceptor.bind(endpoint);
      acceptor.listen();

      _portNumber = acceptor.local_endpoint().port();
    }

    std::string address(ipAddress + ":" + std::to_string(_portNumber));
    std::string addressFileName("." + nameRequester + "-" + nameAcceptor +
                                ".address");

    Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);

    ScopedPublisher p(addressFileName);

    p.write(address);

    // std::cout << address << std::endl;

    preciceDebug("Accept connection at " << address);

    PtrSocket socket(new tcp::socket(*_ioService));

    acceptor.accept(*socket);

    _isConnected = true;

    int remoteSize = requesterCommunicatorSize;

    preciceCheck(remoteSize > 0,
                 "acceptConnection()",
                 "Requester communicator "
                     << "size has to be > 0!");

    _remoteCommunicatorSize = remoteSize;

    _sockets.resize(_remoteCommunicatorSize);
    _clientQueryBuffers.resize(_remoteCommunicatorSize);

    int remoteRank = 0;

    _sockets[remoteRank] = socket;

    send(remoteRank, remoteRank);

    send(0, remoteRank);
    send(1, remoteRank);

    for (int i = 1; i < _remoteCommunicatorSize; ++i) {
      remoteRank = i;

      socket = PtrSocket(new Socket(*_ioService));

      acceptor.accept(*socket); // Waits until connection

      preciceCheck(remoteSize == _remoteCommunicatorSize,
                   "acceptConnection()",
                   "Remote communicator sizes are inconsistent!");
      preciceCheck(_sockets[remoteRank].use_count() == 0,
                   "acceptConnection()",
                   "Duplicate request to connect by same rank (" << remoteRank
                                                                 << ")!");
      _sockets[remoteRank] = socket;

      send(remoteRank, remoteRank);

      send(0, remoteRank);
      send(1, remoteRank);
    }

    acceptor.close();
  } catch (std::exception& e) {
    preciceError("acceptConnection()",
                 "Accepting connection at port " << _portNumber
                                                 << " failed: " << e.what());
  }

  _queryWork = PtrWork(new asio::io_service::work(*_ioService));
  _queryThread = std::thread(&SocketCommunication::onThreadRun, this);
}

void
SocketCommunication::requestConnection(std::string const& nameAcceptor,
                                       std::string const& nameRequester,
                                       int requesterProcessRank,
                                       int requesterCommunicatorSize) {
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);

  using asio::ip::tcp;

  try {
    std::string address;
    std::string addressFileName("." + nameRequester + "-" + nameAcceptor +
                                ".address");

    Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);

    Publisher p(addressFileName);

    p.read(address);

    // std::cout << address << std::endl;

    preciceDebug("Request connection to " << address);

    std::string ipAddress = address.substr(0, address.find(":"));
    std::string portNumber = address.substr(
        ipAddress.length() + 1, address.length() - ipAddress.length() - 1);

    _portNumber = static_cast<unsigned short>(std::stoi(portNumber));

    PtrSocket socket(new Socket(*_ioService));
    tcp::resolver::query query(tcp::v4(), ipAddress, portNumber);

    while (not isConnected()) {
      tcp::resolver resolver(*_ioService);

      {
        tcp::resolver::endpoint_type endpoint = *(resolver.resolve(query));
        boost::system::error_code error = asio::error::host_not_found;
        socket->connect(endpoint, error);

        _isConnected = not error;
      }

      if (not isConnected()) {
        // Wait a little, since after a couple of ten-thousand trials the system
        // seems to get confused and the requester connects wrongly to itself.
        boost::asio::deadline_timer timer(*_ioService,
                                          boost::posix_time::milliseconds(100));
        timer.wait();
      }
    }

    _sockets.push_back(socket);

    send(requesterProcessRank, 0);
    send(requesterCommunicatorSize, 0);

    // Activates sending of queries before actual content. Has to be done after
    // calls to send and before calls to receive.
    _isClient = true;

    _rank = requesterProcessRank;

    int remoteSize = 0;
    int remoteRank = -1;

    receive(remoteRank, 0);
    receive(remoteSize, 0);

    preciceCheck(remoteRank == 0,
                 "requestConnection()",
                 "Acceptor base rank "
                     << "has to be 0 but is " << remoteRank << "!");
    preciceCheck(remoteSize == 1,
                 "requestConnection()",
                 "Acceptor communicator "
                     << "size has to be == 1!");

    _remoteCommunicatorSize = remoteSize;
  } catch (std::exception& e) {
    preciceError("requestConnection()",
                 "Requesting connection failed: " << e.what());
  }

  // NOTE:
  // Keep IO service running so that it fires asynchronous handlers from another
  // thread.
  _queryWork = PtrWork(new asio::io_service::work(*_ioService));
  _queryThread = std::thread([this]() { _ioService->run(); });
}

int
SocketCommunication::requestConnectionAsClient(
    std::string const& nameAcceptor, std::string const& nameRequester) {
  preciceTrace2("requestConnectionAsClient()", nameAcceptor, nameRequester);

  using asio::ip::tcp;

  try {
    std::string address;
    std::string addressFileName("." + nameRequester + "-" + nameAcceptor +
                                ".address");

    Publisher::ScopedChangePrefixDirectory scpd(_addressDirectory);

    Publisher p(addressFileName);

    p.read(address);

    // std::cout << address << std::endl;

    preciceDebug("Request connection to " << address);

    std::string ipAddress = address.substr(0, address.find(":"));
    std::string portNumber = address.substr(
        ipAddress.length() + 1, address.length() - ipAddress.length() - 1);

    _portNumber = static_cast<unsigned short>(std::stoi(portNumber));

    PtrSocket socket(new Socket(*_ioService));
    tcp::resolver::query query(tcp::v4(), ipAddress, portNumber);

    while (not isConnected()) {
      tcp::resolver resolver(*_ioService);

      {
        tcp::resolver::endpoint_type endpoint = *(resolver.resolve(query));
        boost::system::error_code error = asio::error::host_not_found;
        socket->connect(endpoint, error);

        _isConnected = not error;
      }

      if (not isConnected()) {
        // Wait a little, since after a couple of ten-thousand trials the system
        // seems to get confused and the requester connects wrongly to itself.
        boost::asio::deadline_timer timer(*_ioService,
                                          boost::posix_time::milliseconds(100));
        timer.wait();
      }
    }

    _sockets.push_back(socket);

    _isClient = true;

    receive(_rank, 0);

    int remoteSize = 0;
    int remoteRank = -1;

    receive(remoteRank, 0);
    receive(remoteSize, 0);

    preciceCheck(remoteRank == 0,
                 "requestConnection()",
                 "Acceptor base rank "
                     << "has to be 0 but is " << remoteRank << "!");
    preciceCheck(remoteSize == 1,
                 "requestConnection()",
                 "Acceptor communicator "
                     << "size has to be == 1!");

    _remoteCommunicatorSize = remoteSize;
  } catch (std::exception& e) {
    preciceError("requestConnection()",
                 "Requesting connection failed: " << e.what());
  }

  // NOTE:
  // Keep IO service running so that it fires asynchronous handlers from another
  // thread.
  _queryWork = PtrWork(new asio::io_service::work(*_ioService));
  _queryThread = std::thread([this]() { _ioService->run(); });

  return _rank;
}

void
SocketCommunication::closeConnection() {
  preciceTrace("closeConnection()");

  if (not isConnected())
    return;

  if (_queryThread.joinable()) {
    preciceDebug("Shutting down query work, io service, and thread");
    _queryWork.reset();
    _ioService->stop();
    _queryThread.join();
  }

  for (PtrSocket& socket : _sockets) {
    assertion(socket->is_open());
    socket->shutdown(Socket::shutdown_both);
    socket->close();
  }

  _isConnected = false;
}

void
SocketCommunication::startSendPackage(int rankReceiver) {
}

void
SocketCommunication::finishSendPackage() {
}

int
SocketCommunication::startReceivePackage(int rankSender) {
  preciceTrace1("startReceivePackage()", rankSender);
  return rankSender;
}

void
SocketCommunication::finishReceivePackage() {
}

void
SocketCommunication::send(std::string const& itemToSend, int rankReceiver) {
  preciceTrace2("send(string)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver,
             _sockets.size());
  assertion(isConnected());
  size_t size = itemToSend.size() + 1;
  try {
    sendQuery(rankReceiver);
    asio::write(*_sockets[rankReceiver],
                asio::buffer((void*)&size, sizeof(size_t)));
    asio::write(*_sockets[rankReceiver],
                asio::buffer(itemToSend.c_str(), size));
  } catch (std::exception& e) {
    preciceError("send(string)", "Send failed: " << e.what());
  }
}

void
SocketCommunication::send(int* itemsToSend, int size, int rankReceiver) {
  preciceTrace2("send(int*)", size, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver,
             _sockets.size());
  assertion(isConnected());
  try {
    sendQuery(rankReceiver);
    asio::write(*_sockets[rankReceiver],
                asio::buffer((void*)itemsToSend, size * sizeof(int)));
  } catch (std::exception& e) {
    preciceError("send(int*)", "Send failed: " << e.what());
  }
}

Request::SharedPointer
SocketCommunication::aSend(int* itemsToSend, int size, int rankReceiver) {
  preciceTrace2("aSend(int*)", size, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver,
             _sockets.size());
  assertion(isConnected());

  Request::SharedPointer request(new SocketRequest);

  try {
    sendQuery(rankReceiver);
    asio::async_write(*_sockets[rankReceiver],
                      asio::buffer((void*)itemsToSend, size * sizeof(int)),
                      [request](boost::system::error_code const&, std::size_t) {
      static_cast<SocketRequest*>(request.get())->complete();
    });
  } catch (std::exception& e) {
    preciceError("aSend(int*)", "Send failed: " << e.what());
  }

  return request;
}

void
SocketCommunication::send(double* itemsToSend, int size, int rankReceiver) {
  preciceTrace2("send(double*)", size, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver,
             _sockets.size());
  assertion(isConnected());
  try {
    sendQuery(rankReceiver);
    asio::write(*_sockets[rankReceiver],
                asio::buffer((void*)itemsToSend, size * sizeof(double)));
  } catch (std::exception& e) {
    preciceError("send(double*)", "Send failed: " << e.what());
  }
}

Request::SharedPointer
SocketCommunication::aSend(double* itemsToSend, int size, int rankReceiver) {
  preciceTrace2("aSend(double*)", size, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver,
             _sockets.size());
  assertion(isConnected());

  Request::SharedPointer request(new SocketRequest);

  try {
    sendQuery(rankReceiver);
    asio::async_write(*_sockets[rankReceiver],
                      asio::buffer((void*)itemsToSend, size * sizeof(double)),
                      [request](boost::system::error_code const&, std::size_t) {
      static_cast<SocketRequest*>(request.get())->complete();
    });
  } catch (std::exception& e) {
    preciceError("aSend(double*)", "Send failed: " << e.what());
  }

  return request;
}

void
SocketCommunication::send(double itemToSend, int rankReceiver) {
  preciceTrace2("send(double)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver,
             _sockets.size());
  assertion(isConnected());
  try {
    sendQuery(rankReceiver);
    asio::write(*_sockets[rankReceiver],
                asio::buffer((void*)&itemToSend, sizeof(double)));
  } catch (std::exception& e) {
    preciceError("send(double)", "Send failed: " << e.what());
  }
}

Request::SharedPointer
SocketCommunication::aSend(double* itemToSend, int rankReceiver) {
  return aSend(itemToSend, 1, rankReceiver);
}

void
SocketCommunication::send(int itemToSend, int rankReceiver) {
  preciceTrace2("send(int)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver,
             _sockets.size());
  assertion(isConnected());
  try {
    sendQuery(rankReceiver);
    asio::write(*_sockets[rankReceiver],
                asio::buffer((void*)&itemToSend, sizeof(int)));
  } catch (std::exception& e) {
    preciceError("send(int)", "Send failed: " << e.what());
  }
}

Request::SharedPointer
SocketCommunication::aSend(int* itemToSend, int rankReceiver) {
  return aSend(itemToSend, 1, rankReceiver);
}

void
SocketCommunication::send(bool itemToSend, int rankReceiver) {
  preciceTrace2("send(bool)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver,
             _sockets.size());
  assertion(isConnected());
  try {
    sendQuery(rankReceiver);
    asio::write(*_sockets[rankReceiver],
                asio::buffer((void*)&itemToSend, sizeof(bool)));
  } catch (std::exception& e) {
    preciceError("send(double)", "Send failed: " << e.what());
  }
}

Request::SharedPointer
SocketCommunication::aSend(bool* itemToSend, int rankReceiver) {
  preciceTrace1("aSend(bool*)", rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver,
             _sockets.size());
  assertion(isConnected());

  Request::SharedPointer request(new SocketRequest);

  try {
    sendQuery(rankReceiver);
    asio::async_write(*_sockets[rankReceiver],
                      asio::buffer((void*)itemToSend, sizeof(bool)),
                      [request](boost::system::error_code const&, std::size_t) {
      static_cast<SocketRequest*>(request.get())->complete();
    });
  } catch (std::exception& e) {
    preciceError("aSend(bool*)", "Send failed: " << e.what());
  }

  return request;
}

int
SocketCommunication::receive(std::string& itemToReceive, int rankSender) {
  preciceTrace1("receive(string)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2((rankSender >= 0) && (rankSender < (int)_sockets.size()),
             rankSender,
             _sockets.size());
  assertion(isConnected());
  size_t size = 0;
  try {
    asio::read(*_sockets[rankSender],
               asio::buffer((void*)&size, sizeof(size_t)));
    char* msg = new char[size];
    asio::read(*_sockets[rankSender], asio::buffer((void*)msg, size));
    itemToReceive = msg;
    delete[] msg;
  } catch (std::exception& e) {
    preciceError("receive(string)", "Receive failed: " << e.what());
  }
  preciceDebug("Received " << itemToReceive << " from rank " << rankSender);
  if (isServer())
    receiveNextQuery(rankSender);
  return rankSender;
}

int
SocketCommunication::receive(int* itemsToReceive, int size, int rankSender) {
  preciceTrace2("receive(int*)", size, rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2((rankSender >= 0) && (rankSender < (int)_sockets.size()),
             rankSender,
             _sockets.size());
  assertion(isConnected());
  try {
    asio::read(*_sockets[rankSender],
               asio::buffer((void*)itemsToReceive, size * sizeof(int)));
  } catch (std::exception& e) {
    preciceError("receive(int*)", "Receive failed: " << e.what());
  }
  preciceDebug("Received from rank " << rankSender);
  if (isServer())
    receiveNextQuery(rankSender);
  return rankSender;
}

int
SocketCommunication::receive(double* itemsToReceive, int size, int rankSender) {
  preciceTrace2("receive(double*)", size, rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2((rankSender >= 0) && (rankSender < (int)_sockets.size()),
             rankSender,
             _sockets.size());
  assertion(isConnected());
  try {
    asio::read(*_sockets[rankSender],
               asio::buffer((void*)itemsToReceive, size * sizeof(double)));
  } catch (std::exception& e) {
    preciceError("receive(double*)", "Receive failed: " << e.what());
  }
  preciceDebug("Received from rank " << rankSender);
  if (isServer())
    receiveNextQuery(rankSender);
  return rankSender;
}

int
SocketCommunication::receive(double& itemToReceive, int rankSender) {
  preciceTrace1("receive(double)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2((rankSender >= 0) && (rankSender < (int)_sockets.size()),
             rankSender,
             _sockets.size());
  assertion(isConnected());
  try {
    asio::read(*_sockets[rankSender],
               asio::buffer((void*)&itemToReceive, sizeof(double)));
  } catch (std::exception& e) {
    preciceError("receive(double)", "Receive failed: " << e.what());
  }
  preciceDebug("Received " << itemToReceive << " from rank " << rankSender);
  if (isServer())
    receiveNextQuery(rankSender);
  return rankSender;
}

int
SocketCommunication::receive(int& itemToReceive, int rankSender) {
  preciceTrace1("receive(int)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2((rankSender >= 0) && (rankSender < (int)_sockets.size()),
             rankSender,
             _sockets.size());
  assertion(isConnected());
  try {
    asio::read(*_sockets[rankSender],
               asio::buffer((void*)&itemToReceive, sizeof(int)));
  } catch (std::exception& e) {
    preciceError("receive(int)", "Receive failed: " << e.what());
  }
  preciceDebug("Received " << itemToReceive << " from rank " << rankSender);
  if (isServer())
    receiveNextQuery(rankSender);
  return rankSender;
}

int
SocketCommunication::receive(bool& itemToReceive, int rankSender) {
  preciceTrace1("receive(bool)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2((rankSender >= 0) && (rankSender < (int)_sockets.size()),
             rankSender,
             _sockets.size());
  assertion(isConnected());
  try {
    asio::read(*_sockets[rankSender],
               asio::buffer((void*)&itemToReceive, sizeof(bool)));
  } catch (std::exception& e) {
    preciceError("receive(bool)", "Receive failed: " << e.what());
  }
  preciceDebug("Received " << itemToReceive << " from rank " << rankSender);
  if (isServer())
    receiveNextQuery(rankSender);
  return rankSender;
}

int
SocketCommunication::getSenderRank(int desiredRank) {
  preciceTrace1("getSenderRank()", desiredRank);
  int chosenRank = -1;
  if (isClient()) {
    preciceDebug("Client process uses server process 0");
    assertion1(desiredRank != ANY_SENDER, desiredRank);
    chosenRank = desiredRank;
  }
  while (chosenRank == -1) {
    preciceDebug("Server process looks in client queries");
    std::unique_lock<std::mutex> lock(_requestMutex);
    if (not _clientQueries.empty()) {
      if (desiredRank == ANY_SENDER) {
        chosenRank = *_clientQueries.begin();
        _clientQueries.erase(_clientQueries.begin());
      } else {
        std::set<int>::iterator iter = _clientQueries.find(desiredRank);
        if (iter != _clientQueries.end()) {
          chosenRank = *iter;
          _clientQueries.erase(iter);
        }
      }
    }
    if (chosenRank == -1) {
      _requestCondition.wait(lock);
    }
  }
  preciceDebug("Chose sender rank to be " << chosenRank);
  return chosenRank;
}

void
SocketCommunication::sendQuery(int receiverRank) {
  preciceTrace1("sendQuery()", receiverRank);
  if (isClient()) {
    preciceDebug("Process is a client, send query");
    asio::write(*_sockets[receiverRank],
                asio::buffer((void*)&_rank, sizeof(int)));
  }
}

void
SocketCommunication::receiveNextQuery(int senderRank) {
  preciceTrace1("receiveNextQuery()", senderRank);
  assertion2((senderRank >= 0) && (senderRank < _remoteCommunicatorSize),
             senderRank,
             _remoteCommunicatorSize);
  PtrSocket& socket = _sockets[senderRank];
  int* buffer = &_clientQueryBuffers[senderRank];
  int bufferSize = sizeof(int);
  asio::async_read(
      *socket,
      asio::buffer((void*)buffer, bufferSize),
      boost::bind(
          &SocketCommunication::onAsyncReceive, this, _1, _2, senderRank));
}

void
SocketCommunication::onThreadRun() {
  preciceTrace("onThreadRun()");
  for (size_t i = 0; i < _sockets.size(); i++) {
    receiveNextQuery(i);
  }
  _ioService->run(); // Enable asynchronous callbacks
}

void
SocketCommunication::onAsyncReceive(const boost::system::error_code& error,
                                    size_t bytesTransferred,
                                    int clientIndex) {
  preciceTrace2("onAsyncReceive()", bytesTransferred, clientIndex);
  if (error) {
    preciceDebug("Aborting receiving request from client " << clientIndex);
    return;
  }
  assertion2((clientIndex >= 0) && (clientIndex < _remoteCommunicatorSize),
             clientIndex,
             _remoteCommunicatorSize);
  assertion2(_clientQueryBuffers[clientIndex] == clientIndex,
             _clientQueryBuffers[clientIndex],
             clientIndex);
  std::unique_lock<std::mutex> lock(_requestMutex);
  _clientQueries.insert(clientIndex);
  _requestCondition.notify_one();
}

bool
SocketCommunication::isClient() {
  return _isClient;
}

bool
SocketCommunication::isServer() {
  return !_isClient;
}

std::string
SocketCommunication::getIpAddress() {
  preciceTrace("getIpAddress()");
  std::ostringstream oss;

#ifdef _WIN32
  oss << "127.0.0.1";
#else
  int querySocket = socket(AF_INET, SOCK_STREAM, 0);

  assertion(querySocket >= 0);

  preciceDebug("Looking for IP address of network \"" << _networkName << "\"");

  struct ifreq request;
  struct if_nameindex* nameInterface = if_nameindex();

  assertion(nameInterface);

  struct if_nameindex* itNameInterface = nameInterface;

  while (itNameInterface && itNameInterface->if_name) {
    preciceCheck(strlen(itNameInterface->if_name) < IFNAMSIZ,
                 "acceptConnection()",
                 "Network interface \" " << itNameInterface->if_name
                                         << "\" has too long name");

    strncpy(request.ifr_name,
            itNameInterface->if_name,
            IFNAMSIZ); // Copy interface name

    int ifNameLength = strlen(itNameInterface->if_name);

    request.ifr_name[ifNameLength] = 0; // Add C-string 0

    if (ioctl(querySocket, SIOCGIFADDR, &request) >= 0) {
      preciceDebug(itNameInterface->if_name
                   << ": " << inet_ntoa(((struct sockaddr_in*)&request.ifr_addr)
                                            ->sin_addr));

      if (strcmp(itNameInterface->if_name, _networkName.c_str()) == 0) {
        oss << inet_ntoa(((struct sockaddr_in*)&request.ifr_addr)->sin_addr);
      }
    } else {
      preciceCheck(strcmp(itNameInterface->if_name, _networkName.c_str()) != 0,
                   "acceptConnection()",
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
}
} // namespace precice, com

#endif // not PRECICE_NO_SOCKETS
