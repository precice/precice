#ifndef PRECICE_NO_SOCKETS

#pragma once

#include "com/Communication.hpp"
#include <boost/asio.hpp>
#include "logging/Logger.hpp"
#include <thread>
#include "com/SocketSendQueue.hpp"

namespace precice
{
namespace com
{
/// Implements Communication by using sockets.
class SocketCommunication : public Communication
{
public:
  SocketCommunication(unsigned short     portNumber       = 0,
                      bool               reuseAddress     = false,
                      std::string const &networkName      = "lo",
                      std::string const &addressDirectory = ".");

  explicit SocketCommunication(std::string const &addressDirectory);

  virtual ~SocketCommunication();

  virtual size_t getRemoteCommunicatorSize() override;

  virtual void acceptConnection(std::string const &acceptorName,
                                std::string const &requesterName,
                                int                acceptorRank) override;

  virtual void acceptConnectionAsServer(std::string const &acceptorName,
                                        std::string const &requesterName,
                                        int                acceptorRank,
                                        int                requesterCommunicatorSize) override;

  virtual void requestConnection(std::string const &acceptorName,
                                 std::string const &requesterName,
                                 int                requesterRank,
                                 int                requesterCommunicatorSize) override;

  virtual void requestConnectionAsClient(std::string   const &acceptorName,
                                         std::string   const &requesterName,
                                         std::set<int> const &acceptorRanks,
                                         int                  requesterRank) override;


  virtual void closeConnection() override;

  /// Sends a std::string to process with given rank.
  virtual void send(std::string const &itemToSend, int rankReceiver) override;

  /// Sends an array of integer values.
  virtual void send(const int *itemsToSend, int size, int rankReceiver) override;

  /// Asynchronously sends an array of integer values.
  virtual PtrRequest aSend(const int *itemsToSend, int size, int rankReceiver) override;

  /// Sends an array of double values.
  virtual void send(const double *itemsToSend, int size, int rankReceiver) override;

  /// Asynchronously sends an array of double values.
  virtual PtrRequest aSend(const double *itemsToSend, int size, int rankReceiver) override;

  virtual PtrRequest aSend(std::vector<double> const & itemsToSend, int rankReceiver) override;
  
  /// Sends a double to process with given rank.
  virtual void send(double itemToSend, int rankReceiver) override;

  /// Asynchronously sends a double to process with given rank.
  virtual PtrRequest aSend(const double & itemToSend, int rankReceiver) override;

  /// Sends an int to process with given rank.
  virtual void send(int itemToSend, int rankReceiver) override;

  /// Asynchronously sends an int to process with given rank.
  virtual PtrRequest aSend(const int & itemToSend, int rankReceiver) override;

  /// Sends a bool to process with given rank.
  virtual void send(bool itemToSend, int rankReceiver) override;

  /// Asynchronously sends a bool to process with given rank.
  virtual PtrRequest aSend(const bool & itemToSend, int rankReceiver) override;

  /// Receives a std::string from process with given rank.
  virtual void receive(std::string &itemToReceive, int rankSender) override;

  /// Receives an array of integer values.
  virtual void receive(int *itemsToReceive, int size, int rankSender) override;

  /// Receives an array of double values.
  virtual void receive(double *itemsToReceive, int size, int rankSender) override;

  /// Asynchronously receives an array of double values.
  virtual PtrRequest aReceive(double *itemsToReceive,
                              int     size,
                              int     rankSender) override;

  virtual PtrRequest aReceive(std::vector<double> & itemsToReceive, int rankSender) override;
  
  /// Receives a double from process with given rank.
  virtual void receive(double &itemToReceive, int rankSender) override;

  /// Asynchronously receives a double from process with given rank.
  virtual PtrRequest aReceive(double &itemToReceive, int rankSender) override;

  /// Receives an int from process with given rank.
  virtual void receive(int &itemToReceive, int rankSender) override;

  /// Asynchronously receives an int from process with given rank.
  virtual PtrRequest aReceive(int &itemToReceive, int rankSender) override;

  /// Receives a bool from process with given rank.
  virtual void receive(bool &itemToReceive, int rankSender) override;

  /// Asynchronously receives a bool from process with given rank.
  virtual PtrRequest aReceive(bool &itemToReceive, int rankSender) override;

  void send(std::vector<int> const &v, int rankReceiver) override;
  void receive(std::vector<int> &v, int rankSender) override;

  void send(std::vector<double> const &v, int rankReceiver) override;
  void receive(std::vector<double> &v, int rankSender) override;
  
private:
  logging::Logger _log{"com::SocketCommunication"};

  /// Port used for socket connection.
  unsigned short _portNumber;

  bool _reuseAddress;

  /// Name of network to communicate over.
  std::string _networkName;

  /// Directory where IP address is exchanged by file.
  std::string _addressDirectory;

  using IOService     = boost::asio::io_service;
  using TCP           = boost::asio::ip::tcp;
  using SocketService = boost::asio::stream_socket_service<TCP>;
  using Socket        = boost::asio::basic_stream_socket<TCP, SocketService>;
  using Work          = boost::asio::io_service::work;
  
  std::shared_ptr<IOService> _ioService;
  std::shared_ptr<Work> _work;
  std::thread _thread;
  
  /// Remote rank -> socket map
  std::map<int, std::shared_ptr<Socket>> _sockets;

  SendQueue _queue;

  bool isClient();
  bool isServer();

  std::string getIpAddress();
};
} // namespace com
} // namespace precice

#endif // not PRECICE_NO_SOCKETS
