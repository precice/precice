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

  ~SocketCommunication() override;

  size_t getRemoteCommunicatorSize() override;

  void acceptConnection(std::string const &acceptorName,
                        std::string const &requesterName,
                        std::string const &tag,
                        int                acceptorRank,
                        int                rankOffset = 0) override;

  void acceptConnectionAsServer(std::string const &acceptorName,
                                std::string const &requesterName,
                                std::string const &tag,
                                int                acceptorRank,
                                int                requesterCommunicatorSize) override;

  void requestConnection(std::string const &acceptorName,
                         std::string const &requesterName,
                         std::string const &tag,
                         int                requesterRank,
                         int                requesterCommunicatorSize) override;

  void requestConnectionAsClient(std::string const &  acceptorName,
                                 std::string const &  requesterName,
                                 std::string const &  tag,
                                 std::set<int> const &acceptorRanks,
                                 int                  requesterRank) override;

  void closeConnection() override;

  /// Sends a std::string to process with given rank.
  void send(std::string const &itemToSend, Rank rankReceiver) override;

  /// Sends an array of integer values.
  void send(precice::span<const int> itemsToSend, Rank rankReceiver) override;

  /// Asynchronously sends an array of integer values.
  PtrRequest aSend(precice::span<const int> itemsToSend, Rank rankReceiver) override;

  /// Sends an array of double values.
  void send(precice::span<const double> itemsToSend, Rank rankReceiver) override;

  /// Asynchronously sends an array of double values.
  PtrRequest aSend(precice::span<const double> itemsToSend, Rank rankReceiver) override;

  /// Sends a double to process with given rank.
  void send(double itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends a double to process with given rank.
  PtrRequest aSend(const double &itemToSend, Rank rankReceiver) override;

  /// Sends an int to process with given rank.
  void send(int itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends an int to process with given rank.
  PtrRequest aSend(const int &itemToSend, Rank rankReceiver) override;

  /// Sends a bool to process with given rank.
  void send(bool itemToSend, Rank rankReceiver) override;

  /// Asynchronously sends a bool to process with given rank.
  PtrRequest aSend(const bool &itemToSend, Rank rankReceiver) override;

  /// Receives a std::string from process with given rank.
  void receive(std::string &itemToReceive, Rank rankSender) override;

  /// Receives an array of integer values.
  void receive(precice::span<int> itemsToReceive, Rank rankSender) override;

  /// Receives an array of double values.
  void receive(precice::span<double> itemsToReceive, Rank rankSender) override;

  /// Asynchronously receives an array of double values.
  PtrRequest aReceive(precice::span<double> itemsToReceive,
                      int                   rankSender) override;

  /// Receives a double from process with given rank.
  void receive(double &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives a double from process with given rank.
  PtrRequest aReceive(double &itemToReceive, Rank rankSender) override;

  /// Receives an int from process with given rank.
  void receive(int &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives an int from process with given rank.
  PtrRequest aReceive(int &itemToReceive, Rank rankSender) override;

  /// Receives a bool from process with given rank.
  void receive(bool &itemToReceive, Rank rankSender) override;

  /// Asynchronously receives a bool from process with given rank.
  PtrRequest aReceive(bool &itemToReceive, Rank rankSender) override;

  void prepareEstablishment(std::string const &acceptorName,
                            std::string const &requesterName) override;

  void cleanupEstablishment(std::string const &acceptorName,
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

  using IOContext = boost::asio::io_context;
  using Socket    = boost::asio::ip::tcp::socket;
  using WorkGuard = boost::asio::executor_work_guard<IOContext::executor_type>;

  std::shared_ptr<IOContext> _ioContext;
  std::unique_ptr<WorkGuard> _workGuard;
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
