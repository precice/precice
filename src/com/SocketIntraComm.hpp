#pragma once

#include <memory>
#include <string>

#include "com/IntraCommunication.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"

namespace precice::com {

class SocketCommunication;

/**
 * @brief Socket-based implementation of IntraCommunication.
 *
 * Uses TCP sockets for all operations. Collective operations (broadcast,
 * reduce, allreduce) are built from point-to-point socket communication.
 *
 * Internally composes a SocketCommunication for the actual socket I/O,
 * but does NOT inherit from Communication. The public interface is purely
 * IntraCommunication.
 */
class SocketIntraComm : public IntraCommunication {
public:
  SocketIntraComm(unsigned short portNumber       = 0,
                  bool           reuseAddress     = false,
                  std::string    networkName      = {},
                  std::string    addressDirectory = ".");

  ~SocketIntraComm() override;

  /// @name Connection
  /// @{
  bool isConnected() override;
  size_t getRemoteCommunicatorSize() override;
  void connectIntraComm(std::string const &participantName,
                        std::string const &tag,
                        int                rank,
                        int                size) override;
  void closeConnection() override;
  /// @}

  /// @name Reduction (built from point-to-point)
  /// @{
  void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive) override;
  void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank) override;
  void reduceSum(int itemToSend, int &itemToReceive) override;
  void reduceSum(int itemToSend, int &itemToReceive, Rank primaryRank) override;

  void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive) override;
  void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank) override;
  void allreduceSum(double itemToSend, double &itemToReceive) override;
  void allreduceSum(double itemToSend, double &itemToReceive, Rank primaryRank) override;
  void allreduceSum(int itemToSend, int &itemToReceive) override;
  void allreduceSum(int itemToSend, int &itemToReceive, Rank primaryRank) override;
  /// @}

  /// @name Broadcast (built from point-to-point)
  /// @{
  void broadcast(precice::span<const int> itemsToSend) override;
  void broadcast(precice::span<int> itemsToReceive, Rank rankBroadcaster) override;
  void broadcast(int itemToSend) override;
  void broadcast(int &itemToReceive, Rank rankBroadcaster) override;
  void broadcast(precice::span<const double> itemsToSend) override;
  void broadcast(precice::span<double> itemsToReceive, Rank rankBroadcaster) override;
  void broadcast(double itemToSend) override;
  void broadcast(double &itemToReceive, Rank rankBroadcaster) override;
  void broadcast(bool itemToSend) override;
  void broadcast(bool &itemToReceive, Rank rankBroadcaster) override;
  void broadcast(std::vector<int> const &v) override;
  void broadcast(std::vector<int> &v, Rank rankBroadcaster) override;
  void broadcast(std::vector<double> const &v) override;
  void broadcast(std::vector<double> &v, Rank rankBroadcaster) override;
  /// @}

  /// @name Point-to-point (socket-based)
  /// @{
  void send(std::string const &itemToSend, Rank rankReceiver) override;
  void send(precice::span<const int> itemsToSend, Rank rankReceiver) override;
  PtrRequest aSend(precice::span<const int> itemsToSend, Rank rankReceiver) override;
  void send(precice::span<const double> itemsToSend, Rank rankReceiver) override;
  PtrRequest aSend(precice::span<const double> itemsToSend, Rank rankReceiver) override;
  void send(double itemToSend, Rank rankReceiver) override;
  PtrRequest aSend(const double &itemToSend, Rank rankReceiver) override;
  void send(int itemToSend, Rank rankReceiver) override;
  PtrRequest aSend(const int &itemToSend, Rank rankReceiver) override;
  void send(bool itemToSend, Rank rankReceiver) override;
  PtrRequest aSend(const bool &itemToSend, Rank rankReceiver) override;

  void receive(std::string &itemToReceive, Rank rankSender) override;
  void receive(precice::span<int> itemsToReceive, Rank rankSender) override;
  void receive(precice::span<double> itemsToReceive, Rank rankSender) override;
  PtrRequest aReceive(precice::span<double> itemsToReceive, int rankSender) override;
  void receive(double &itemToReceive, Rank rankSender) override;
  PtrRequest aReceive(double &itemToReceive, Rank rankSender) override;
  void receive(int &itemToReceive, Rank rankSender) override;
  PtrRequest aReceive(int &itemToReceive, Rank rankSender) override;
  void receive(bool &itemToReceive, Rank rankSender) override;
  PtrRequest aReceive(bool &itemToReceive, Rank rankSender) override;
  /// @}

private:
  logging::Logger _log{"com::SocketIntraComm"};

  /// Internal socket communication used for the actual I/O.
  std::unique_ptr<SocketCommunication> _socketComm;

  /// Number of secondary ranks
  size_t _remoteCommunicatorSize = 0;

  /// The rank offset used by the internal SocketCommunication
  static constexpr int RANK_OFFSET = 1;

  /// Helper to iterate over remote communicator ranks
  size_t remoteSize() const { return _remoteCommunicatorSize; }
};

} // namespace precice::com
