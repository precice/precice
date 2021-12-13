#ifndef PRECICE_NO_MPI

#pragma once

#include <mpi.h>
#include <set>
#include <stddef.h>
#include <string>

#include "MPICommunication.hpp"
#include "logging/Logger.hpp"
#include "precice/types.hpp"
#include "utils/Parallel.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace com {
/**
 * @brief Provides connection methods for processes located in one communicator.
 *
 * This communication class can be used when the communicating participants are
 * either compiled into one executable or, are started by one mpi execution call
 * on the command line.
 *
 * It is imporant, that all processes in the used communicator have to
 * participate in the communication. If one of the processes does not call
 * either acceptConnection(), or closeConnection(), a deadlock is achieved.
 */
class MPIDirectCommunication : public MPICommunication {
public:
  /** Creates a MPI Direct communication based on the current global communicator.
   */
  MPIDirectCommunication();

  virtual ~MPIDirectCommunication();

  /**
   * @brief Returns the number of processes in the remote communicator.
   *
   * @pre A connection to the remote participant has been setup.
   */
  virtual size_t getRemoteCommunicatorSize() override;

  /** See precice::com::Communication::acceptConnection().
   * @attention Calls precice::utils::Parallel::splitCommunicator()
   * if local and global communicators are equal.
   */
  virtual void acceptConnection(std::string const &acceptorName,
                                std::string const &requesterName,
                                std::string const &tag,
                                int                acceptorRank,
                                int                rankOffset = 0) override;

  virtual void acceptConnectionAsServer(std::string const &acceptorName,
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
  virtual void requestConnection(std::string const &acceptorName,
                                 std::string const &requesterName,
                                 std::string const &tag,
                                 int                requesterRank,
                                 int                requesterCommunicatorSize) override;

  virtual void requestConnectionAsClient(std::string const &  acceptorName,
                                         std::string const &  requesterName,
                                         std::string const &  tag,
                                         std::set<int> const &acceptorRanks,
                                         int                  requesterRank) override
  {
    PRECICE_ASSERT(false, "Not implemented!");
  }

  /// See precice::com::Communication::closeConnection().
  virtual void closeConnection() override;

  virtual void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank rankMaster) override;

  virtual void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive) override;

  virtual void reduceSum(int itemToSend, int &itemsToReceive, Rank rankMaster) override;

  virtual void reduceSum(int itemToSend, int &itemsToReceive) override;

  virtual void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank rankMaster) override;

  virtual void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive) override;

  virtual void allreduceSum(double itemToSend, double &itemsToReceive, Rank rankMaster) override;

  virtual void allreduceSum(double itemToSend, double &itemsToReceive) override;

  virtual void allreduceSum(int itemToSend, int &itemsToReceive, Rank rankMaster) override;

  virtual void allreduceSum(int itemToSend, int &itemsToReceive) override;

  virtual void broadcast(precice::span<const int> itemsToSend) override;

  virtual void broadcast(precice::span<int> itemsToReceive, Rank rankBroadcaster) override;

  virtual void broadcast(int itemToSend) override;

  virtual void broadcast(int &itemToReceive, Rank rankBroadcaster) override;

  virtual void broadcast(precice::span<const double> itemsToSend) override;

  virtual void broadcast(precice::span<double> itemsToReceive, Rank rankBroadcaster) override;

  virtual void broadcast(double itemToSend) override;

  virtual void broadcast(double &itemToReceive, Rank rankBroadcaster) override;

  virtual void broadcast(bool itemToSend) override;

  virtual void broadcast(bool &itemToReceive, Rank rankBroadcaster) override;

private:
  virtual MPI_Comm &communicator(Rank rank = 0) override;

  virtual Rank rank(int rank) override;

  logging::Logger _log{"com::MPIDirectCommunication"};

  /// CommState to use
  utils::Parallel::CommStatePtr _commState = nullptr;

protected:
  /// Turn the rank adjustment into a noop for direct communication
  virtual int adjustRank(Rank rank) const override;
};

} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
