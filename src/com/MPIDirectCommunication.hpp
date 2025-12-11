#ifndef PRECICE_NO_MPI

#pragma once

#include <mpi.h>
#include <set>
#include <stddef.h>
#include <string>

#include "com/MPICommunication.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"
#include "utils/Parallel.hpp"
#include "utils/assertion.hpp"

namespace precice::com {
/**
 * @brief Provides connection methods for processes located in one communicator.
 *
 * This communication class can be used when the communicating participants are
 * either compiled into one executable or, are started by one mpi execution call
 * on the command line.
 *
 * It is important, that all processes in the used communicator have to
 * participate in the communication. If one of the processes does not call
 * either acceptConnection(), or closeConnection(), a deadlock is achieved.
 */
class MPIDirectCommunication : public MPICommunication {
public:
  /** Creates a MPI Direct communication based on the current global communicator.
   */
  MPIDirectCommunication();

  ~MPIDirectCommunication() override;

  /**
   * @brief Returns the number of processes in the remote communicator.
   *
   * @pre A connection to the remote participant has been setup.
   */
  size_t getRemoteCommunicatorSize() override;

  /** See precice::com::Communication::acceptConnection().
   * @attention Calls precice::utils::Parallel::splitCommunicator()
   * if local and global communicators are equal.
   */
  void acceptConnection(std::string const &acceptorName,
                        std::string const &requesterName,
                        std::string const &tag,
                        int                acceptorRank,
                        int                rankOffset = 0) override;

  void acceptConnectionAsServer(std::string const &acceptorName,
                                std::string const &requesterName,
                                std::string const &tag,
                                int                acceptorRank,
                                int                requesterCommunicatorSize) override
  {
    PRECICE_ASSERT(false, "Not implemented!");
  }

  /** See precice::com::Communication::requestConnection().
   * @attention Calls precice::utils::Parallel::splitCommunicator()
   * if local and global communicators are equal.
   */
  void requestConnection(std::string const &acceptorName,
                         std::string const &requesterName,
                         std::string const &tag,
                         int                requesterRank,
                         int                requesterCommunicatorSize) override;

  void requestConnectionAsClient(std::string const   &acceptorName,
                                 std::string const   &requesterName,
                                 std::string const   &tag,
                                 std::set<int> const &acceptorRanks,
                                 int                  requesterRank) override
  {
    PRECICE_ASSERT(false, "Not implemented!");
  }

  /// See precice::com::Communication::closeConnection().
  void closeConnection() override;

  void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank) override;

  void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive) override;

  void reduceSum(int itemToSend, int &itemsToReceive, Rank primaryRank) override;

  void reduceSum(int itemToSend, int &itemsToReceive) override;

  void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank primaryRank) override;

  void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive) override;

  void allreduceSum(double itemToSend, double &itemsToReceive, Rank primaryRank) override;

  void allreduceSum(double itemToSend, double &itemsToReceive) override;

  void allreduceSum(int itemToSend, int &itemsToReceive, Rank primaryRank) override;

  void allreduceSum(int itemToSend, int &itemsToReceive) override;

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

private:
  MPI_Comm &communicator(Rank rank = 0) override;

  Rank rank(int rank) override;

  logging::Logger _log{"com::MPIDirectCommunication"};

  /// CommState to use
  utils::Parallel::CommStatePtr _commState = nullptr;

protected:
  /// Turn the rank adjustment into a noop for direct communication
  int adjustRank(Rank rank) const override;
};

} // namespace precice::com

#endif // not PRECICE_NO_MPI
