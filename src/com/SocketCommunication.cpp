// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
# ifndef PRECICE_NO_SOCKETS

#include "SocketCommunication.hpp"
#include <sstream>
#include <fstream>
#include <boost/asio.hpp>
#include <boost/bind.hpp>

namespace precice {
namespace com {

namespace asio = boost::asio;

tarch::logging::Log SocketCommunication:: _log("precice::com::SocketCommunication");

SocketCommunication:: SocketCommunication
(
  const std::string& network,
  int                port,
  const std::string& ipExchangeDirectory )
:
  _network(network),
  _port(port),
  _ipExchangeDirectory(ipExchangeDirectory),
  _processRank(-1),
  _isConnected(false),
  _remoteCommunicatorSize(0),
  _ioService(new IOService()),
  _sockets(),
  _queryWork(),
  _queryThread(),
  _clientQueries(),
  _clientQueryBuffers()
{
  _rankOffset = 0;
}

SocketCommunication:: ~SocketCommunication()
{
  if (_isConnected){
    closeConnection();
  }
}

bool SocketCommunication:: isConnected()
{
  return _isConnected;
}

int SocketCommunication:: getRemoteCommunicatorSize()
{
  assertion(_isConnected);
  return _remoteCommunicatorSize;
}

void SocketCommunication:: acceptConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                acceptorProcessRank,
  int                acceptorCommunicatorSize )
{
  preciceTrace2("acceptConnection()", nameAcceptor, nameRequester);
  preciceCheck ( acceptorCommunicatorSize == 1, "acceptConnection()",
                 "Acceptor of socket connection can only have one process!" );
  std::string ipFilename(_ipExchangeDirectory + "." + nameRequester + "-portname");
  using asio::ip::tcp;
  try {
    std::ostringstream address;
    // Query for IP address
    int querySocket = socket(AF_INET, SOCK_STREAM, 0);
    assertion(querySocket >= 0);
    preciceDebug("Looking for IP address of network \"" << _network << "\"");
    struct ifreq request;
    struct if_nameindex* nameInterface = if_nameindex();
    assertion(nameInterface);
    struct if_nameindex* itNameInterface = nameInterface;
    while (itNameInterface && itNameInterface->if_name){
      preciceCheck(strlen(itNameInterface->if_name) < IFNAMSIZ,
                   "acceptConnection()", "Network interface \" "
                   << itNameInterface->if_name << "\" has too long name");
      strncpy(request.ifr_name, itNameInterface->if_name, IFNAMSIZ); // Copy interface name
      int ifNameLength = strlen(itNameInterface->if_name);
      request.ifr_name[ifNameLength] = 0; // Add C-string 0
      if (ioctl(querySocket, SIOCGIFADDR, &request) >= 0){
        preciceDebug(itNameInterface->if_name << ": "
                     << inet_ntoa(((struct sockaddr_in*) &request.ifr_addr)->sin_addr));
        if(strcmp(itNameInterface->if_name, _network.c_str()) == 0){
          address << inet_ntoa(((struct sockaddr_in*) &request.ifr_addr)->sin_addr);
        }
      }
      else {
        preciceCheck(strcmp(itNameInterface->if_name, _network.c_str()) != 0,
                     "acceptConnection()", "Could not obtain network IP from "
                     << "network \"" << itNameInterface->if_name << "\"");
      }
      itNameInterface++;
    }
    if_freenameindex(nameInterface);
    close(querySocket);
    preciceCheck(not address.str().empty(), "acceptConnection()",
                 "Network \"" << _network << "\" not found for socket connection!");

    // Write server address to file
    preciceDebug("Writing server ip address to file " << ipFilename);
    std::ofstream outFile;
    outFile.open ((ipFilename +  "~").c_str(), std::ios::out);
    outFile << address.str();
    outFile.close();
    // To give the file first a different name prevents early reading errors
    rename( (ipFilename + "~").c_str(), ipFilename.c_str() );

    preciceDebug("Accept connection at " << address.str() << ":" << _port);
    tcp::acceptor acceptor(*_ioService, tcp::endpoint(tcp::v4(), _port));
    PtrSocket socket ( new tcp::socket(*_ioService) );
    acceptor.accept(*socket); // Waits until connection
    _isConnected = true;
    int remoteRank = -1;
    int remoteSize = 0;
    asio::read ( *socket, asio::buffer((void*)&remoteRank, sizeof(int)) );
    asio::read ( *socket, asio::buffer((void*)&remoteSize, sizeof(int)) );
    preciceCheck ( remoteSize > 0, "acceptConnection()", "Requester communicator "
                   << "size has to be > 0!" );
    _remoteCommunicatorSize = remoteSize;
    preciceDebug("Received rank=" << remoteRank << ", size=" << remoteSize);
    _sockets.resize(_remoteCommunicatorSize);
    _clientQueryBuffers.resize(_remoteCommunicatorSize);
    _sockets[remoteRank] = socket;
    send ( acceptorProcessRank, remoteRank );
    send ( acceptorCommunicatorSize, remoteRank );
    for ( int i=1; i < _remoteCommunicatorSize; i++ ){ // Connect to remaining processes
      socket = PtrSocket(new Socket(*_ioService));
      acceptor.accept(*socket); // Waits until connection
      asio::read ( *socket, asio::buffer((void*)&remoteRank, sizeof(int)) );
      asio::read ( *socket, asio::buffer((void*)&remoteSize, sizeof(int)) );
      preciceCheck(remoteSize == _remoteCommunicatorSize, "acceptConnection()",
                   "Remote communicator sizes are inconsistent!" );
      preciceCheck(_sockets[remoteRank].use_count() == 0, "acceptConnection()",
                   "Duplicate request to connect by same rank (" << remoteRank << ")!");
      _sockets[remoteRank] = socket;
      send ( acceptorProcessRank, remoteRank );
      send ( acceptorCommunicatorSize, remoteRank );
    }
    acceptor.close();
  }
  catch (std::exception& e){
    preciceError("acceptConnection()", "Accepting connection at port " << _port
                 << " failed: " << e.what());
  }

  //if ( utils::Parallel::getLocalProcessRank() == 0 ){
  if (remove(ipFilename.c_str()) != 0){
    preciceWarning("acceptConnection()", "Could not remove ip address file \""
                   << ipFilename << "\"!" );
  }

  _queryWork = PtrWork(new asio::io_service::work(*_ioService));
  _queryThread = boost::thread(&SocketCommunication::onThreadRun, this);
}

void SocketCommunication:: requestConnection
(
  const std::string& nameAcceptor,
  const std::string& nameRequester,
  int                requesterProcessRank,
  int                requesterCommunicatorSize )
{
  preciceTrace2("requestConnection()", nameAcceptor, nameRequester);
  using asio::ip::tcp;
  try {
    // Read server address from file
    std::string ipFilename(_ipExchangeDirectory + "." + nameRequester + "-portname");
    preciceDebug("Reading server ip address from file " << ipFilename);
    std::ifstream inFile;
    do {
      inFile.open(ipFilename.c_str(), std::ios::in);
    } while (not inFile);
    std::string serverAddress;
    inFile >> serverAddress;
    inFile.close();
    preciceDebug("Read connection info \"" << serverAddress
                 << "\" from file " << ipFilename);

    PtrSocket socket(new Socket(*_ioService));
    std::ostringstream portStream;
    portStream << _port;
    tcp::resolver::query query(tcp::v4(), serverAddress.c_str(), portStream.str().c_str());
    while (not _isConnected){ // since resolver does not wait until server is up
      tcp::resolver resolver(*_ioService);
      tcp::resolver::endpoint_type endpoint = *(resolver.resolve(query));
      preciceDebug("Request connection at " << endpoint);
      boost::system::error_code error = asio::error::host_not_found;
      socket->connect(endpoint, error);
      _isConnected = not error;
      if (not _isConnected){
        // Wait a little, since after a couple of ten-thousand trials the system
        // seems to get confused and the requester connects wrongly to itself.
        boost::asio::deadline_timer timer(*_ioService, boost::posix_time::milliseconds(100));
        timer.wait();
      }
    }
    _sockets.push_back(socket); // Only one socket in requrestConnection possible

    send(requesterProcessRank, 0);
    send(requesterCommunicatorSize, 0);
    int remoteSize = 0;
    int remoteRank = -1;
    // Activates sending of queries before actual content. Has to be done after
    // calls to send.
    _processRank = requesterProcessRank;
    receive(remoteRank, 0);
    receive(remoteSize, 0);
    preciceCheck(remoteRank == 0, "requestConnection()", "Acceptor base rank "
                 << "has to be 0 but is " << remoteRank << "!");
    preciceCheck(remoteSize == 1, "requestConnection()", "Acceptor communicator "
                 << "size has to be == 1!");
    _remoteCommunicatorSize = remoteSize;
  }
  catch (std::exception& e){
    preciceError("requestConnection()", "Requesting connection failed: " << e.what());
  }
}

void SocketCommunication:: closeConnection()
{
  preciceTrace("closeConnection()");
  assertion(_isConnected);
  if (_queryThread.joinable()){
    preciceDebug("Shutting down query work, io service, and thread");
    _queryWork.reset();
    _ioService->stop();
    _queryThread.join();
  }
  for (PtrSocket& socket : _sockets ) {
    assertion(socket->is_open());
    socket->shutdown(Socket::shutdown_both);
    socket->close();
  }
  _isConnected = false;
}

void SocketCommunication:: startSendPackage
(
  int rankReceiver )
{

}

void SocketCommunication:: finishSendPackage()
{

}

int SocketCommunication:: startReceivePackage
(
  int rankSender )
{
  preciceTrace1("startReceivePackage()", rankSender);
  return rankSender;
}

void SocketCommunication:: finishReceivePackage()
{

}

void SocketCommunication:: send
(
  const std::string& itemToSend,
  int                rankReceiver )
{
  preciceTrace2("send(string)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2 ( (rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
               rankReceiver, _sockets.size() );
  assertion(_isConnected);
  size_t size = itemToSend.size() + 1;
  try {
    sendQuery(rankReceiver);
    asio::write ( *_sockets[rankReceiver], asio::buffer((void*)&size, sizeof(size_t)) );
    asio::write ( *_sockets[rankReceiver], asio::buffer(itemToSend.c_str(), size) );
  }
  catch (std::exception& e){
    preciceError("send(string)", "Send failed: " << e.what());
  }
}

void SocketCommunication:: send
(
  int* itemsToSend,
  int  size,
  int  rankReceiver )
{
  preciceTrace2("send(int*)", size, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2((rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
             rankReceiver, _sockets.size());
  assertion(_isConnected);
  try {
    sendQuery(rankReceiver);
    asio::write(*_sockets[rankReceiver],
                asio::buffer((void*)itemsToSend, size*sizeof(int)));
  }
  catch (std::exception& e){
    preciceError("send(int*)", "Send failed: " << e.what());
  }
}

void SocketCommunication:: send
(
  double* itemsToSend,
  int     size,
  int     rankReceiver )
{
  preciceTrace2("send(double*)", size, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2 ( (rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
               rankReceiver, _sockets.size() );
  assertion(_isConnected);
  try {
    sendQuery(rankReceiver);
    asio::write(*_sockets[rankReceiver],
                asio::buffer((void*)itemsToSend, size*sizeof(double)));
  }
  catch (std::exception& e){
    preciceError("send(double*)", "Send failed: " << e.what());
  }
}

void SocketCommunication:: send
(
  double itemToSend,
  int    rankReceiver )
{
  preciceTrace2("send(double)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2 ( (rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
               rankReceiver, _sockets.size() );
  assertion(_isConnected);
  try {
    sendQuery(rankReceiver);
    asio::write ( *_sockets[rankReceiver],
                  asio::buffer((void*)&itemToSend, sizeof(double)) );
  }
  catch (std::exception& e){
    preciceError("send(double)", "Send failed: " << e.what());
  }
}

void SocketCommunication:: send
(
  int itemToSend,
  int rankReceiver )
{
  preciceTrace2("send(int)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2 ( (rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
               rankReceiver, _sockets.size() );
  assertion(_isConnected);
  try {
    sendQuery(rankReceiver);
    asio::write ( *_sockets[rankReceiver],
                  asio::buffer((void*)&itemToSend, sizeof(int)) );
  }
  catch (std::exception& e){
    preciceError("send(double)", "Send failed: " << e.what());
  }
}

void SocketCommunication:: send
(
  bool itemToSend,
  int  rankReceiver )
{
  preciceTrace2("send(bool)", itemToSend, rankReceiver);
  rankReceiver = rankReceiver - _rankOffset;
  assertion2 ( (rankReceiver >= 0) && (rankReceiver < (int)_sockets.size()),
               rankReceiver, _sockets.size() );
  assertion(_isConnected);
  try {
    sendQuery(rankReceiver);
    asio::write ( *_sockets[rankReceiver],
                  asio::buffer((void*)&itemToSend, sizeof(bool)) );
  }
  catch (std::exception& e){
    preciceError("send(double)", "Send failed: " << e.what());
  }
}

int SocketCommunication:: receive
(
  std::string& itemToReceive,
  int          rankSender )
{
  preciceTrace1("receive(string)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2 ( (rankSender >= 0) && (rankSender < (int)_sockets.size()),
               rankSender, _sockets.size() );
  assertion(_isConnected);
  size_t size = 0;
  try {
    asio::read ( *_sockets[rankSender], asio::buffer((void*)&size, sizeof(size_t)) );
    char* msg = new char[size];
    asio::read ( *_sockets[rankSender], asio::buffer((void*)msg, size) );
    itemToReceive = msg;
    delete[] msg;
  }
  catch (std::exception& e){
    preciceError("receive(string)", "Receive failed: " << e.what());
  }
  preciceDebug("Received " << itemToReceive << " from rank " << rankSender);
  if (_processRank == -1) receiveNextQuery(rankSender);
  return rankSender;
}

int SocketCommunication:: receive
(
  int* itemsToReceive,
  int  size,
  int  rankSender )
{
  preciceTrace2("receive(int*)", size, rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2 ( (rankSender >= 0) && (rankSender < (int)_sockets.size()),
               rankSender, _sockets.size() );
  assertion(_isConnected);
  try {
    asio::read ( *_sockets[rankSender],
                 asio::buffer((void*)itemsToReceive, size*sizeof(int)) );
  }
  catch (std::exception& e){
    preciceError("receive(int*)", "Receive failed: " << e.what());
  }
  preciceDebug("Received from rank " << rankSender);
  if (_processRank == -1) receiveNextQuery(rankSender);
  return rankSender;
}

int SocketCommunication:: receive
(
  double* itemsToReceive,
  int     size,
  int     rankSender )
{
  preciceTrace2("receive(double*)", size, rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2 ( (rankSender >= 0) && (rankSender < (int)_sockets.size()),
               rankSender, _sockets.size() );
  assertion(_isConnected);
  try {
    asio::read ( *_sockets[rankSender],
                 asio::buffer((void*)itemsToReceive, size*sizeof(double)) );
  }
  catch (std::exception& e){
    preciceError("receive(double*)", "Receive failed: " << e.what());
  }
  preciceDebug("Received from rank " << rankSender);
  if (_processRank == -1) receiveNextQuery(rankSender);
  return rankSender;
}

int SocketCommunication:: receive
(
  double& itemToReceive,
  int     rankSender )
{
  preciceTrace1("receive(double)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2 ( (rankSender >= 0) && (rankSender < (int)_sockets.size()),
               rankSender, _sockets.size() );
  assertion(_isConnected);
  try {
    asio::read ( *_sockets[rankSender],
                 asio::buffer((void*)&itemToReceive, sizeof(double)) );
  }
  catch (std::exception& e){
    preciceError("receive(double)", "Receive failed: " << e.what());
  }
  preciceDebug("Received " << itemToReceive << " from rank " << rankSender);
  if (_processRank == -1) receiveNextQuery(rankSender);
  return rankSender;
}

int SocketCommunication:: receive
(
  int& itemToReceive,
  int  rankSender )
{
  preciceTrace1("receive(int)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2 ( (rankSender >= 0) && (rankSender < (int)_sockets.size()),
               rankSender, _sockets.size() );
  assertion(_isConnected);
  try {
    asio::read ( *_sockets[rankSender],
                 asio::buffer((void*)&itemToReceive, sizeof(int)) );
  }
  catch (std::exception& e){
    preciceError("receive(int)", "Receive failed: " << e.what());
  }
  preciceDebug("Received " << itemToReceive << " from rank " << rankSender);
  if (_processRank == -1) receiveNextQuery(rankSender);
  return rankSender;
}

int SocketCommunication:: receive
(
  bool& itemToReceive,
  int   rankSender )
{
  preciceTrace1("receive(bool)", rankSender);
  rankSender = rankSender - _rankOffset;
  rankSender = getSenderRank(rankSender);
  assertion2 ( (rankSender >= 0) && (rankSender < (int)_sockets.size()),
               rankSender, _sockets.size() );
  assertion(_isConnected);
  try {
    asio::read ( *_sockets[rankSender],
                 asio::buffer((void*)&itemToReceive, sizeof(bool)) );
  }
  catch (std::exception& e){
    preciceError("receive(bool)", "Receive failed: " << e.what());
  }
  preciceDebug("Received " << itemToReceive << " from rank " << rankSender);
  if (_processRank == -1) receiveNextQuery(rankSender);
  return rankSender;
}

int SocketCommunication:: getSenderRank
(
  int desiredRank )
{
  preciceTrace1("getSenderRank()", desiredRank);
  int chosenRank = -1;
  if (_processRank != -1){
    preciceDebug("Client process uses server process 0");
    assertion1(desiredRank != ANY_SENDER, desiredRank);
    chosenRank = desiredRank;
  }
  while (chosenRank == -1){
    preciceDebug("Server process looks in client queries");
    boost::mutex::scoped_lock lock(_requestMutex);
    if (not _clientQueries.empty()){
      if (desiredRank == ANY_SENDER){
        chosenRank = *_clientQueries.begin();
        _clientQueries.erase(_clientQueries.begin());
      }
      else {
        std::set<int>::iterator iter = _clientQueries.find(desiredRank);
        if (iter != _clientQueries.end()){
          chosenRank = *iter;
          _clientQueries.erase(iter);
        }
      }
    }
    if (chosenRank == -1){
      _requestCondition.wait(lock);
    }
  }
  preciceDebug("Chose sender rank to be " << chosenRank);
  return chosenRank;
}

void SocketCommunication:: sendQuery
(
  int receiverRank )
{
  preciceTrace1("sendQuery()", receiverRank);
  if (_processRank != -1){
    preciceDebug("Process is a client, send query");
    asio::write ( *_sockets[receiverRank],
                  asio::buffer((void*)&_processRank, sizeof(int)) );
  }
}

void SocketCommunication:: receiveNextQuery
(
  int senderRank )
{
  preciceTrace1("receiveNextQuery()", senderRank);
  assertion2((senderRank >= 0) && (senderRank < _remoteCommunicatorSize),
             senderRank, _remoteCommunicatorSize);
  PtrSocket& socket = _sockets[senderRank];
  int* buffer = &_clientQueryBuffers[senderRank];
  int bufferSize = sizeof(int);
  asio::async_read(*socket, asio::buffer((void*)buffer, bufferSize),
      boost::bind(&SocketCommunication::onAsyncReceive, this, _1, _2, senderRank));
}

void SocketCommunication:: onThreadRun()
{
  preciceTrace("onThreadRun()");
  for (size_t i=0; i < _sockets.size(); i++){
    receiveNextQuery(i);
  }
  _ioService->run(); // Enable asynchronous callbacks
}

void SocketCommunication:: onAsyncReceive
(
  const boost::system::error_code& error,
  size_t                           bytesTransferred,
  int                              clientIndex )
{
  preciceTrace2("onAsyncReceive()", bytesTransferred, clientIndex);
  if (error){
    preciceDebug("Aborting receiving request from client " << clientIndex);
    return;
  }
  assertion2((clientIndex >= 0) && (clientIndex < _remoteCommunicatorSize),
             clientIndex, _remoteCommunicatorSize);
  assertion2(_clientQueryBuffers[clientIndex] == clientIndex,
             _clientQueryBuffers[clientIndex], clientIndex);
  boost::mutex::scoped_lock lock(_requestMutex);
  _clientQueries.insert(clientIndex);
  _requestCondition.notify_one();
}

}} // namespace precice, com

#endif // not PRECICE_NO_SOCKETS
