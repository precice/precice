// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#ifndef PRECICE_COM_MPI_DIRECT_COMMUNICATION_HPP_
#define PRECICE_COM_MPI_DIRECT_COMMUNICATION_HPP_

#include "MPICommunication.hpp"

#include "tarch/logging/Log.h"

#include <string>
#include <vector>

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
  /**
   * @brief Constructor.
   */
  MPIDirectCommunication();

  /**
   * @brief Destructor.
   */
  virtual ~MPIDirectCommunication();

  /**
   * @brief Returns true, if a connection to a remote participant has been
   * setup.
   */
  virtual bool
  isConnected() {
    return _isConnected;
  }

  /**
   * @brief Returns the number of processes in the remote communicator.
   *
   * Precondition: a connection to the remote participant has been setup.
   */
  virtual int getRemoteCommunicatorSize();

  /**
   * @brief See precice::com::Communication::acceptConnection().
   */
  virtual void acceptConnection(std::string const& nameAcceptor,
                                std::string const& nameRequester,
                                int acceptorProcessRank,
                                int acceptorCommunicatorSize);

  virtual void
  acceptConnectionAsServer(std::string const& nameAcceptor,
                           std::string const& nameRequester,
                           int requesterCommunicatorSize) {
    preciceError("acceptConnectionAsServer()", "Not implemented!");
  }

  /**
   * @brief See precice::com::Communication::requestConnection().
   */
  virtual void requestConnection(std::string const& nameAcceptor,
                                 std::string const& nameRequester,
                                 int requesterProcessRank,
                                 int requesterCommunicatorSize);

  virtual int
  requestConnectionAsClient(std::string const& nameAcceptor,
                            std::string const& nameRequester) {
    preciceError("requestConnectionAsClient()", "Not implemented!");
  }

  /**
   * @brief See precice::com::Communication::closeConnection().
   */
  virtual void closeConnection();

  virtual void broadcast();

  virtual void broadcast(int* itemsToSend, int size);

  virtual void broadcast(int* itemsToReceive, int size, int rankBroadcaster);

  virtual void broadcast(int itemToSend);

  virtual void broadcast(int& itemToReceive, int rankBroadcaster);

  virtual void broadcast(double* itemsToSend, int size);

  virtual void broadcast(double* itemsToReceive, int size, int rankBroadcaster);

  virtual void broadcast(double itemToSend);

  virtual void broadcast(double& itemToReceive, int rankBroadcaster);

  virtual void broadcast(bool itemToSend);

  virtual void broadcast(bool& itemToReceive, int rankBroadcaster);

private:
  virtual MPI_Comm& communicator(int rank = 0);

  virtual int rank(int rank);

  // @brief Logging device.
  static tarch::logging::Log _log;

  MPI_Comm _communicator;

  // @brief Global communicator, as given by utils::Parallel::getDefaultComm().
  MPI_Comm _globalCommunicator;

  // @brief Communicator for communicator between process groups.
  MPI_Comm _localCommunicator;

  bool _isConnected;

  /**
   * @brief Returns ID belonging to a group of processes.
   *
   * Erroneous, if called before exchangeGroupInformation.
   */
  int getGroupID(std::string const& accessorName);

  /**
   * @brief Returns rank of leading process of a group.
   *
   * Erroneous, if called before exchangeGroupInformation.
   */
  int getLeaderRank(std::string const& accessorName);
};
}
} // namespace precice, com

#endif /* PRECICE_COM_MPI_DIRECT_COMMUNICATION_HPP_ */

#endif // not PRECICE_NO_MPI
