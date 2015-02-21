// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_M2N_POINT_TO_POINT_COMMUNICATION_HPP_
#define PRECICE_M2N_POINT_TO_POINT_COMMUNICATION_HPP_

#include "DistributedCommunication.hpp"

#include "com/Communication.hpp"
#include "com/CommunicationFactory.hpp"
#include "tarch/logging/Log.h"
#include "com/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace m2n {
/**
 * @brief Point-to-point communication implementation of
 *        DistributedCommunication.
 *
 * Direct communication of local data subsets is performed between processes of
 * coupled participants. The two supported implementations of direct
 * communication are SocketCommunication and MPIPortsCommunication which can be
 * supplied via their corresponding instantiation factories
 * SocketCommunicationFactory and MPIPortsCommunicationFactory.
 *
 * For the detailed implementation documentation refer to
 * PointToPointCommunication.cpp.
 */
class PointToPointCommunication : public DistributedCommunication {
public:
  /**
   * @brief Constructor.
   */
  PointToPointCommunication(com::PtrCommunicationFactory communicationFactory,
                            mesh::PtrMesh mesh);

  /**
   * @brief Destructor.
   */
  virtual ~PointToPointCommunication();

  /**
   * @brief Returns true, if a connection to a remote participant has been
   *        established.
   */
  virtual bool isConnected();

  /**
   * @brief Accepts connection from participant, which has to call
   *        requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  virtual void acceptConnection(std::string const& nameAcceptor,
                                std::string const& nameRequester);

  /**
   * @brief Requests connection from participant, which has to call
   *        acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  virtual void requestConnection(std::string const& nameAcceptor,
                                 std::string const& nameRequester);

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection();

  /**
   * @brief Sends a subset of local double values corresponding to local indices
   *        deduced from the current and remote vertex distributions.
   */
  virtual void send(double* itemsToSend, int size, int valueDimension = 1);

  /**
   * @brief Receives a subset of local double values corresponding to local
   *        indices deduced from the current and remote vertex distributions.
   */
  virtual void receive(double* itemsToReceive,
                       int size,
                       int valueDimension = 1);

private:
  static tarch::logging::Log _log;

  com::PtrCommunicationFactory _communicationFactory;

  std::map<int, com::PtrCommunication> _communications;

  /**
   * @brief Local (for process rank in the current participant) communication
   *        map that defines a mapping from a process rank in the remote
   *        participant to an array of local data indices, which define a subset
   *        of local (for process rank in the current participant) data to be
   *        communicated between the current process rank and the remote process
   *        rank.
   *
   * Example. Assume that the current process rank is 3. Assume that its
   * `_communicationMap' is
   *
   *   1 -> {1, 3}
   *   4 -> {0, 2}
   *
   * then it means that the current process (with rank 3)
   * - has to communicate (send/receive) data with local indices 1 and 3 with
   *   the remote process with rank 1;
   * - has to communicate (send/receive) data with local indices 0 and 2 with
   *   the remote process with rank 4.
   */
  std::map<int, std::vector<int>> _communicationMap;

  bool _isConnected;
  bool _isAcceptor;
};
}
} // namespace precice, m2n

#endif /* PRECICE_M2N_POINT_TO_POINT_COMMUNICATION_HPP_ */
