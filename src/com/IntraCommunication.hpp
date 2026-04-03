#pragma once

#include <stddef.h>
#include <string>
#include <vector>

#include "com/Request.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"
#include "precice/span.hpp"

namespace precice::com {

/**
 * @brief Interface for intra-participant communication (primary-secondary).
 *
 * This is a standalone interface that replaces the use of the generic
 * Communication class for intra-participant communication. It provides:
 * - Efficient collective operations (broadcast, reduce, allreduce)
 * - Point-to-point send/receive for partition data exchange
 * - Simplified connection setup for primary-secondary topology
 *
 * Rank convention: rank 0 is the primary, ranks 1..N-1 are secondaries.
 * No rank offset adjustment is needed.
 *
 * @see MPIIntraComm for an MPI-based implementation using native collectives.
 * @see SocketIntraComm for a socket-based implementation.
 */
class IntraCommunication {
public:
  IntraCommunication &operator=(IntraCommunication &&) = delete;

  virtual ~IntraCommunication() = default;

  /// @name Connection Setup
  /// @{

  /// Returns true if a connection has been established.
  virtual bool isConnected() = 0;

  /**
   * @brief Returns the number of secondary ranks in the intra-communicator.
   *
   * For MPI, this is the communicator size minus 1 (excluding primary).
   * For sockets, this is the number of connected secondaries.
   *
   * @pre A connection has been set up via connectIntraComm().
   */
  virtual size_t getRemoteCommunicatorSize() = 0;

  /**
   * @brief Establishes the intra-participant communication.
   *
   * Sets up communication between the primary rank (rank 0) and all
   * secondary ranks (1..size-1).
   *
   * @param[in] participantName Name of the participant.
   * @param[in] tag Tag for establishing this connection.
   * @param[in] rank The current rank in the participant (0 = primary).
   * @param[in] size Total number of ranks in the participant.
   */
  virtual void connectIntraComm(std::string const &participantName,
                                std::string const &tag,
                                int                rank,
                                int                size) = 0;

  /// Closes the intra-participant communication.
  virtual void closeConnection() = 0;

  /// @}

  /// @name Reduction
  /// @{

  /// Performs a reduce summation on the primary rank. Called by the primary.
  virtual void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive) = 0;

  /// Performs a reduce summation, sending to the given primary rank. Called by secondaries.
  virtual void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank) = 0;

  /// Performs a reduce summation of a single int on the primary rank.
  virtual void reduceSum(int itemToSend, int &itemToReceive) = 0;

  /// Performs a reduce summation of a single int, sending to the given primary rank.
  virtual void reduceSum(int itemToSend, int &itemToReceive, Rank primaryRank) = 0;

  /// Performs an allreduce summation on the primary rank, then distributes result.
  virtual void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive) = 0;

  /// Performs an allreduce summation, for secondary ranks.
  virtual void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank) = 0;

  /// Performs an allreduce summation of a single double on the primary rank.
  virtual void allreduceSum(double itemToSend, double &itemToReceive) = 0;

  /// Performs an allreduce summation of a single double, for secondary ranks.
  virtual void allreduceSum(double itemToSend, double &itemToReceive, Rank primaryRank) = 0;

  /// Performs an allreduce summation of a single int on the primary rank.
  virtual void allreduceSum(int itemToSend, int &itemToReceive) = 0;

  /// Performs an allreduce summation of a single int, for secondary ranks.
  virtual void allreduceSum(int itemToSend, int &itemToReceive, Rank primaryRank) = 0;

  /// @}

  /// @name Broadcast
  /// @{

  virtual void broadcast(precice::span<const int> itemsToSend) = 0;
  virtual void broadcast(precice::span<int> itemsToReceive, Rank rankBroadcaster) = 0;

  virtual void broadcast(int itemToSend) = 0;
  virtual void broadcast(int &itemToReceive, Rank rankBroadcaster) = 0;

  virtual void broadcast(precice::span<const double> itemsToSend) = 0;
  virtual void broadcast(precice::span<double> itemsToReceive, Rank rankBroadcaster) = 0;

  virtual void broadcast(double itemToSend) = 0;
  virtual void broadcast(double &itemToReceive, Rank rankBroadcaster) = 0;

  virtual void broadcast(bool itemToSend) = 0;
  virtual void broadcast(bool &itemToReceive, Rank rankBroadcaster) = 0;

  virtual void broadcast(std::vector<int> const &v) = 0;
  virtual void broadcast(std::vector<int> &v, Rank rankBroadcaster) = 0;

  virtual void broadcast(std::vector<double> const &v) = 0;
  virtual void broadcast(std::vector<double> &v, Rank rankBroadcaster) = 0;

  /// @}

  /// @name Send
  /// @{

  virtual void send(std::string const &itemToSend, Rank rankReceiver) = 0;
  virtual void send(precice::span<const int> itemsToSend, Rank rankReceiver) = 0;
  virtual PtrRequest aSend(precice::span<const int> itemsToSend, Rank rankReceiver) = 0;
  virtual void send(precice::span<const double> itemsToSend, Rank rankReceiver) = 0;
  virtual PtrRequest aSend(precice::span<const double> itemsToSend, Rank rankReceiver) = 0;
  virtual void send(double itemToSend, Rank rankReceiver) = 0;
  virtual PtrRequest aSend(const double &itemToSend, Rank rankReceiver) = 0;
  virtual void send(int itemToSend, Rank rankReceiver) = 0;
  virtual PtrRequest aSend(const int &itemToSend, Rank rankReceiver) = 0;
  virtual void send(bool itemToSend, Rank rankReceiver) = 0;
  virtual PtrRequest aSend(const bool &itemToSend, Rank rankReceiver) = 0;

  /// @}

  /// @name Receive
  /// @{

  virtual void receive(std::string &itemToReceive, Rank rankSender) = 0;
  virtual void receive(precice::span<int> itemsToReceive, Rank rankSender) = 0;
  virtual void receive(precice::span<double> itemsToReceive, Rank rankSender) = 0;
  virtual PtrRequest aReceive(precice::span<double> itemsToReceive, int rankSender) = 0;
  virtual void receive(double &itemToReceive, Rank rankSender) = 0;
  virtual PtrRequest aReceive(double &itemToReceive, Rank rankSender) = 0;
  virtual void receive(int &itemToReceive, Rank rankSender) = 0;
  virtual PtrRequest aReceive(int &itemToReceive, Rank rankSender) = 0;
  virtual void receive(bool &itemToReceive, Rank rankSender) = 0;
  virtual PtrRequest aReceive(bool &itemToReceive, Rank rankSender) = 0;

  /// @}

  /// @name Range communication
  /// @{

  /** Tag used to specify which type of vector to return.
   * @see receiveRange()
   */
  template <typename T>
  struct AsVectorTag {
  };

  /// Sends a range of doubles (size + content)
  void sendRange(precice::span<const double> itemsToSend, Rank rankReceiver);

  /// Sends a range of ints (size + content)
  void sendRange(precice::span<const int> itemsToSend, Rank rankReceiver);

  /// Receives a range of ints as a vector<int>
  std::vector<int> receiveRange(Rank rankSender, AsVectorTag<int>);

  /// Receives a range of doubles as a vector<double>
  std::vector<double> receiveRange(Rank rankSender, AsVectorTag<double>);

  /// @}

protected:
  bool _isConnected = false;
};

/// Allows to use @ref IntraCommunication::AsVectorTag in a less verbose way.
template <typename T>
inline constexpr auto intraAsVector = IntraCommunication::AsVectorTag<T>{};

} // namespace precice::com
