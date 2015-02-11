// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#ifndef PRECICE_COM_MPI_PORTS_COMMUNICATION_HPP_
#define PRECICE_COM_MPI_PORTS_COMMUNICATION_HPP_

#include "MPICommunication.hpp"

#include "tarch/logging/Log.h"

#include <vector>

namespace precice {
namespace com {

/**
 * @brief Provides connection methods based on MPI ports (part of MPI 2.0).
 *
 * The two participants to be connected can be run in two process groups started
 * up individually, i.e. not within the same process group.
 */
class MPIPortsCommunication : public MPICommunication {
public:
  /**
   * @brief Constructor.
   */
  MPIPortsCommunication(std::string const& addressDirectory = ".");

  /**
   * @brief Destructor.
   */
  virtual ~MPIPortsCommunication();

  /**
   * @brief Returns true, if a connection to a remote participant has been
   * setup.
   */
  virtual bool isConnected();

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

  /**
   * @brief See precice::com::Communication::requestConnection().
   */
  virtual void requestConnection(std::string const& nameAcceptor,
                                 std::string const& nameRequester,
                                 int requesterProcessRank,
                                 int requesterCommunicatorSize);

  /**
   * @brief See precice::com::Communication::closeConnection().
   */
  virtual void closeConnection();

  virtual void send(std::string const& itemToSend, int rankReceiver);

  virtual void send(int* itemsToSend, int size, int rankReceiver);

  virtual void send(double* itemsToSend, int size, int rankReceiver);

  virtual void send(double itemToSend, int rankReceiver);

  virtual void send(int itemToSend, int rankReceiver);

  virtual void send(bool itemToSend, int rankReceiver);

  virtual int receive(std::string& itemToReceive, int rankSender);

  virtual int receive(int* itemsToReceive, int size, int rankSender);

  virtual int receive(double* itemsToReceive, int size, int rankSender);

  virtual int receive(double& itemToReceive, int rankSender);

  virtual int receive(int& itemToReceive, int rankSender);

  virtual int receive(bool& itemToReceive, int rankSender);

private:
  // @brief Logging device.
  static tarch::logging::Log _log;

  std::string _addressDirectory;

  std::vector<MPI_Comm> _communicators;

  // @brief Name of the port used for connection.
  char _port_name[MPI_MAX_PORT_NAME];

  bool _isAcceptor;

  // @brief Flag indicating a connection.
  bool _isConnected;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_MPI_PORTS_COMMUNICATION_HPP_ */

#endif // not PRECICE_NO_MPI
