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
 * @brief Implements GlobalCommunication by using a gathering/scattering
 * methodology.
 * Arrays of data are always gathered and scattered at the master. No direct
 * communication
 * between slaves is used.
 * For more details see m2n/GlobalCommunication.hpp
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
   * setup.
   */
  virtual bool isConnected();

  /**
   * @brief Accepts connection from participant, which has to call
   *        requestConnection().
   *
   * If several connections are going in to a server, the server has to call
   * this method, while the clients have to call requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  virtual void acceptConnection(const std::string& nameAcceptor,
                                const std::string& nameRequester);

  /**
   * @brief Requests connection from participant, which has to call
   *        acceptConnection().
   *
   * If several connections are going in to a server, the clients have to call
   * this method, while the server has to call acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  virtual void requestConnection(const std::string& nameAcceptor,
                                 const std::string& nameRequester);

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection();

  /**
   * @brief Sends an array of double values from all slaves (different for each
   *        slave).
   */
  virtual void send(double* itemsToSend, int size, int valueDimension = 1);

  /**
   * @brief All slaves receive an array of doubles (different for each slave).
   */
  virtual void receive(double* itemsToReceive, int size, int valueDimension = 1);

private:
  static tarch::logging::Log _log;

  com::PtrCommunicationFactory _communicationFactory;

  std::map<int, com::PtrCommunication> _communications;
  std::map<int, std::vector<int>> _senderMap;

  /**
   * @brief global communication is set up or not
   */
  bool _isConnected;
  bool _isAcceptor;
};
}
} // namespace precice, m2n

#endif /* PRECICE_M2N_POINT_TO_POINT_COMMUNICATION_HPP_ */
