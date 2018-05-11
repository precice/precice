#ifndef PRECICE_NO_MPI

#pragma once

#include <string>
#include "MPICommunication.hpp"
#include "logging/Logger.hpp"

namespace precice
{
namespace com
{
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
class MPIDirectCommunication : public MPICommunication
{
public:
  MPIDirectCommunication();

  virtual ~MPIDirectCommunication();

  /**
   * @brief Returns the number of processes in the remote communicator.
   *
   * @pre A connection to the remote participant has been setup.
   */
  virtual size_t getRemoteCommunicatorSize();

  /// See precice::com::Communication::acceptConnection().
  virtual void acceptConnection(std::string const &nameAcceptor,
                                std::string const &nameRequester) override;

  virtual void acceptConnectionAsServer(std::string const &nameAcceptor,
                                        std::string const &nameRequester,
                                        int                requesterCommunicatorSize)
  {
    ERROR("Not implemented!");
  }

  /// See precice::com::Communication::requestConnection().
  virtual void requestConnection(std::string const &nameAcceptor,
                                 std::string const &nameRequester,
                                 int                requesterProcessRank,
                                 int                requesterCommunicatorSize);

  virtual int requestConnectionAsClient(std::string const &nameAcceptor,
                                        std::string const &nameRequester)
  {
    ERROR("Not implemented!");
  }

  /// See precice::com::Communication::closeConnection().
  virtual void closeConnection();

  virtual void reduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster);

  virtual void reduceSum(double *itemsToSend, double *itemsToReceive, int size);

  virtual void reduceSum(int itemToSend, int &itemsToReceive, int rankMaster);

  virtual void reduceSum(int itemToSend, int &itemsToReceive);

  virtual void allreduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster);

  virtual void allreduceSum(double *itemsToSend, double *itemsToReceive, int size);

  virtual void allreduceSum(double itemToSend, double &itemsToReceive, int rankMaster);

  virtual void allreduceSum(double itemToSend, double &itemsToReceive);

  virtual void allreduceSum(int itemToSend, int &itemsToReceive, int rankMaster);

  virtual void allreduceSum(int itemToSend, int &itemsToReceive);

  virtual void broadcast(int *itemsToSend, int size);

  virtual void broadcast(int *itemsToReceive, int size, int rankBroadcaster);

  virtual void broadcast(int itemToSend);

  virtual void broadcast(int &itemToReceive, int rankBroadcaster);

  virtual void broadcast(double *itemsToSend, int size);

  virtual void broadcast(double *itemsToReceive, int size, int rankBroadcaster);

  virtual void broadcast(double itemToSend);

  virtual void broadcast(double &itemToReceive, int rankBroadcaster);

  virtual void broadcast(bool itemToSend);

  virtual void broadcast(bool &itemToReceive, int rankBroadcaster);

private:
  virtual MPI_Comm &communicator(int rank = 0);

  virtual int rank(int rank);

  logging::Logger _log{"com::MPIDirectCommunication"};

  MPI_Comm _communicator;

  /// Global communicator, as given by utils::Parallel::getDefaultComm().
  MPI_Comm _globalCommunicator;

  /// Communicator for communicator between process groups.
  MPI_Comm _localCommunicator;

  /**
   * @brief Returns ID belonging to a group of processes.
   *
   * @pre Call exchangeGroupInformation.
   */
  int getGroupID(std::string const &accessorName);

  /**
   * @brief Returns rank of leading process of a group.
   *
   * @pre Call exchangeGroupInformation.
   */
  int getLeaderRank(std::string const &accessorName);
};

} // namespace com
} // namespace precice

#endif // not PRECICE_NO_MPI
