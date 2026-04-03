#ifndef PRECICE_NO_MPI

#pragma once

#include <mpi.h>
#include <string>

#include "com/IntraCommunication.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"
#include "utils/Parallel.hpp"

namespace precice::com {

/**
 * @brief Flat MPI-based implementation of IntraCommunication.
 *
 * Directly uses an MPI communicator for all operations:
 * - Collectives use MPI_Bcast, MPI_Reduce, MPI_Allreduce
 * - Point-to-point uses MPI_Send/MPI_Recv
 *
 * No intermediate abstraction layers. This is the most efficient
 * backend for intra-participant communication when MPI is available.
 */
class MPIIntraComm : public IntraCommunication {
public:
  MPIIntraComm();

  ~MPIIntraComm() override;

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

  /// @name Reduction (MPI_Reduce / MPI_Allreduce)
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

  /// @name Broadcast (MPI_Bcast)
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

  /// @name Point-to-point (MPI_Send / MPI_Recv)
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
  logging::Logger _log{"com::MPIIntraComm"};

  /// The MPI communicator state, directly used for all operations.
  utils::Parallel::CommStatePtr _commState;
};

} // namespace precice::com

#endif // not PRECICE_NO_MPI
