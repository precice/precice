#pragma once

#include <boost/asio.hpp>
#include <map>
#include <memory>
#include <set>
#include <stddef.h>
#include <string>
#include <thread>

#include "com/Communication.hpp"
#include "com/SharedPointer.hpp"
#include "com/SocketSendQueue.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"
#include "utils/networking.hpp"

namespace precice {
namespace com {
/// Implements Communication by using sockets.
class SocketCommunication : public Communication {
public:
  SocketCommunication(unsigned short portNumber       = 0,
                      bool           reuseAddress     = false,
                      std::string    networkName      = utils::networking::loopbackInterfaceName(),
                      std::string    addressDirectory = ".");

  explicit SocketCommunication(std::string const &addressDirectory);

  virtual ~SocketCommunication();

  virtual size_t getRemoteCommunicatorSize() override;

  virtual void acceptConnection(std::string const &acceptorName,
                                std::string const &requesterName,
                                std::string const &tag,
                                int                acceptorRank,
                                int                rankOffset = 0) override;

  virtual void acceptConnectionAsServer(std::string const &acceptorName,
                                        std::string const &requesterName,
                                        std::string const &tag,
                                        int                acceptorRank,
                                        int                requesterCommunicatorSize) override;

  virtual void requestConnection(std::string const &acceptorName,
                                 std::string const &requesterName,
                                 std::string const &tag,
                                 int                requesterRank,
                                 int                requesterCommunicatorSize) override;

  virtual void requestConnectionAsClient(std::string const &  acceptorName,
                                         std::string const &  requesterName,
                                         std::string const &  tag,
                                         std::set<int> const &acceptorRanks,
                                         int                  requesterRank) override;

  virtual void closeConnection() override;

  /// Sends a std::string to process with given rank.
  virtual void send(std::string const &itemToSend, Rank rankReceiver) override;

  /// Sends an array of integer values.
  virtual void send(precice::span<const int> itemsToSend, Rank rankReceiver) override;

  /// Asynchronously sends an array of integer values.
  virtual PtrRequest aSend(precice::span<const int> itemsToSend, Rank rankReceiver) override;

  /// Sends an array of double values.
  virtual void send(precice::span<const double> itemsToSend, Rank rankReceiver) override;

  /// Asynchronously sends an array of double values.
  virtual PtrRequest aSend(precice::span<const double> itemsToSend, Rank rankReceiver) override;

  /// Sends a double to process with given rank.
  virtual void send(double itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends a double to process with given rank.
  virtual PtrRequest aSend(const double &itemToSend, Rank rankReceiver) override;

  /// Sends an int to process with given rank.
  virtual void send(int itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends an int to process with given rank.
  virtual PtrRequest aSend(const int &itemToSend, Rank rankReceiver) override;

  /// Sends a bool to process with given rank.
  virtual void send(bool itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends a bool to process with given rank.
  virtual PtrRequest aSend(const bool &itemToSend, Rank rankReceiver) override;

  /// Receives a std::string from process with given rank.
  virtual void receive(std::string &itemToReceive, Rank rankSender) override;

  /// Receives an array of integer values.
  virtual void receive(precice::span<int> itemsToReceive, Rank rankSender) override;

  /// Receives an array of double values.
  virtual void receive(precice::span<double> itemsToReceive, Rank rankSender) override;

  /// Asynchronously receives an array of double values.
  virtual PtrRequest aReceive(precice::span<double> itemsToReceive,
                              int                   rankSender) override;

  /// Receives a double from process with given rank.
  virtual void receive(double &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives a double from process with given rank.
  virtual PtrRequest aReceive(double &itemToReceive, Rank rankSender) override;

  /// Receives an int from process with given rank.
  virtual void receive(int &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives an int from process with given rank.
  virtual PtrRequest aReceive(int &itemToReceive, Rank rankSender) override;

  /// Receives a bool from process with given rank.
  virtual void receive(bool &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives a bool from process with given rank.
  virtual PtrRequest aReceive(bool &itemToReceive, Rank rankSender) override;

  virtual void prepareEstablishment(std::string const &acceptorName,
                                    std::string const &requesterName) override;

  virtual void cleanupEstablishment(std::string const &acceptorName,
                                    std::string const &requesterName) override;

private:
  logging::Logger _log{"com::SocketCommunication"};

  /// Port used for socket connection.
  unsigned short _portNumber;

  bool _reuseAddress;

  /// Name of network to communicate over.
  std::string _networkName;

  /// Directory where IP address is exchanged by file.
  std::string _addressDirectory;

  using IOService = boost::asio::io_service;
  using Socket    = boost::asio::ip::tcp::socket;
  using Work      = boost::asio::io_service::work;

  std::shared_ptr<IOService> _ioService;
  std::shared_ptr<Work>      _work;
  std::thread                _thread;

  /// Remote rank -> socket map
  std::map<int, std::shared_ptr<Socket>> _sockets;

  SocketSendQueue _queue;

  bool isClient();
  bool isServer();

  std::string getIpAddress();
};
} // namespace com
} // namespace precice
