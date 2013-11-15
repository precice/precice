// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_SOCKETCOMMUNICATION_HPP_
#define PRECICE_COM_SOCKETCOMMUNICATION_HPP_

#include "com/Communication.hpp"
#include "tarch/logging/Log.h"
#include "utils/PointerVector.hpp"
#include <set>
#include <boost/smart_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/asio/io_service.hpp>

namespace boost {
  namespace asio {
    class io_service;
    namespace ip {
      class tcp;
    }
    template<typename Protocol>
    class stream_socket_service;
    template<typename Protocol, typename StreamSocketService>
    class basic_stream_socket;
  }
  namespace system {
    class error_code;
  }
}


// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace com {

/**
 * @brief Implements Communication by using sockets.
 */
class SocketCommunication : public Communication
{
public:

  /**
   * @brief Constructor.
   */
  SocketCommunication (
    const std::string& network,
    int                port,
    const std::string& addressExchangeDirectory );

  /**
   * @brief Destructor.
   */
  virtual ~SocketCommunication();

  /**
   * @brief Returns true, if a connection to a remote participant has been setup.
   */
  virtual bool isConnected();

  /**
   * @brief Returns the number of processes in the remote communicator.
   *
   * Precondition: a connection to the remote participant has been setup.
   */
  virtual int getRemoteCommunicatorSize();

  /**
   * @brief Connects to another participant, which has to call requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  virtual void acceptConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester,
    int                acceptorProcessRank,
    int                acceptorCommunicatorSize );

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  virtual void requestConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester,
    int                requesterProcessRank,
    int                requesterCommunicatorSize );

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection();

  virtual void startSendPackage ( int rankReceiver );

  virtual void finishSendPackage();

  /**
   * @brief Starts to receive messages from rankSender.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int startReceivePackage ( int rankSender );

  virtual void finishReceivePackage();

  /**
   * @brief Sends a std::string to process with given rank.
   */
  virtual void send (
    const std::string& itemToSend,
    int                rankReceiver );

  /**
   * @brief Sends an array of integer values.
   */
  virtual void send (
    int* itemsToSend,
    int  size,
    int  rankReceiver );

  /**
   * @brief Sends an array of double values.
   */
  virtual void send (
    double* itemsToSend,
    int     size,
    int     rankReceiver );

  /**
   * @brief Sends a double to process with given rank.
   */
  virtual void send (
    double itemToSend,
    int    rankReceiver );

  /**
   * @brief Sends an int to process with given rank.
   */
  virtual void send (
    int itemToSend,
    int rankReceiver );

  /**
   * @brief Sends a bool to process with given rank.
   */
  virtual void send (
    bool itemToSend,
    int  rankReceiver );

  /**
   * @brief Receives a std::string from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive (
    std::string& itemToReceive,
    int          rankSender );

  /**
   * @brief Receives an array of integer values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive (
    int* itemsToReceive,
    int  size,
    int  rankSender );

  /**
   * @brief Receives an array of double values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive (
    double* itemsToReceive,
    int     size,
    int     rankSender );

  /**
   * @brief Receives a double from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive (
    double& itemToReceive,
    int     rankSender );

  /**
   * @brief Receives an int from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive (
    int& itemToReceive,
    int  rankSender );

  /**
   * @brief Receives a bool from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive (
    bool& itemToReceive,
    int   rankSender );

private:

  static tarch::logging::Log _log;

  // @brief Name of network to communicate over.
  std::string _network;

  // @brief Port used for socket connection.
  int _port;

  // @brief Directory where IP address is exchanged by file.
  std::string _ipExchangeDirectory;

  // @brief Local process rank sent as query.
  int _processRank;

  bool _isConnected;

  int _remoteCommunicatorSize;

  typedef boost::asio::io_service IOService;
  boost::shared_ptr<IOService> _ioService;

  typedef boost::asio::ip::tcp TCP;
  typedef boost::asio::stream_socket_service<TCP> SocketService;
  typedef boost::asio::basic_stream_socket<TCP,SocketService> Socket;
  typedef boost::shared_ptr<Socket> PtrSocket;
  std::vector<PtrSocket> _sockets;

  typedef boost::asio::io_service::work Work;
  typedef boost::shared_ptr<Work> PtrWork;
  PtrWork _queryWork;

  // @brief Thread for asynchronously receiving send requests of clients.
  boost::thread _queryThread;

  // @brief Stores socket indices of clients waiting with messages.
  std::set<int> _clientQueries;

  // @brief Buffers for receiving entries in _clientQueries.
  std::vector<int> _clientQueryBuffers;

  boost::mutex _requestMutex;

  boost::condition _requestCondition;

  /**
   * @brief Returns a suitable sender rank to receive from.
   *
   * Uses the _clientQueries to choose a suitable rank from. If no suitable rank
   * is contained, waits until the query thread has received one.
   *
   * If the desiredRank == ANY_RANK, the first rank in _clientQueries is chosen.
   */
  int getSenderRank ( int desiredRank );

  /**
   * @brief If the local process is a client process, sends query to receiverRank.
   *
   * The SocketCommuniation class models communication between a client or
   * several clients and a server process. The roles are defined on setup of
   * communication, the acceptor of the connection is the server, requesters are
   * clients. A client needs to send a request to the server, before he sends
   * the actual data. This method determines whether the process is a client
   * process and sends the query if necesssary.
   */
  void sendQuery ( int receiverRank );

  /**
   * @brief Starts asynchronous receiving of next sender query from senderRank.
   *
   * When a query is received in onAsyncReceive() by a specific sender rank, the
   * next query is not automatically received and has to be triggered by this
   * method.
   */
  void receiveNextQuery ( int senderRank );

  /**
   * @brief Callback on starting to run the query thread.
   */
  void onThreadRun();

  /**
   * @brief Callback on completion of asyncronous receive operation.
   *
   * Is used by query thread to store send queries of clients for master thread.
   */
  void onAsyncReceive (
    const boost::system::error_code& error,
    size_t                           bytesTransferred,
    int                              clientIndex );
};

}} // namespace precice, com

#endif /* PRECICE_COM_SOCKETCOMMUNICATION_HPP_ */
