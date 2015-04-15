// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_M2N_DISTRIBUTED_COMMUNICATION_HPP_
#define PRECICE_M2N_DISTRIBUTED_COMMUNICATION_HPP_

#include "mesh/SharedPointer.hpp"

namespace precice {
namespace m2n {

/**
 * @brief Interface for all distributed solver to solver communication classes.
 *
 *
 * By default, communication is done within the local communication space. In
 * order to connect to a different communication space, i.e. coupling participant,
 * the methods acceptConnection() and requestConnection() have to called by the
 * two participants which intend to establish a connection. All following
 * communication and process ranking refers to the remote communication space
 * afterwards.
 *
 * This interface organizes the communication between 2 distributed participants.
 * The core communication (e.g. Sockets or MPI) is still handled in com/Communication.hpp .
 * The information on how the mesh is distributed can be accessed through the member variable
 * _mesh.
 *
 * The class offers methods to communicate between all processors. This can either
 * mean that data is communicated in a distributed way, in case of arrays, or that single values
 * are broadcasted.
 *
 */
class DistributedCommunication
{
public:
  using SharedPointer = std::shared_ptr<DistributedCommunication>;

public:

  DistributedCommunication(mesh::PtrMesh mesh)
  :
  _mesh(mesh)
  {}

  /**
   * @brief Destructor, empty.
   */
  virtual ~DistributedCommunication() {}

  /**
   * @brief Returns true, if a connection to a remote participant has been setup.
   */
  virtual bool isConnected() =0;

  /**
   * @brief Connects to another participant, which has to call requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  virtual void acceptConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester) =0;

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  virtual void requestConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester) =0;

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection() =0;

  /**
   * @brief Sends an array of double values from all slaves (different for each slave).
   */
  virtual void send (
    double* itemsToSend,
    int     size,
    int     valueDimension) =0;

  /**
   * @brief All slaves receive an array of doubles (different for each slave).
   */
  virtual void receive (
    double* itemsToReceive,
    int     size,
    int     valueDimension) =0;

protected:
  /**
   * @brief mesh that dictates the distribution of this mapping TODO maybe change this directly to vertexDistribution
   */
  mesh::PtrMesh _mesh;
};



}} // namespace precice, m2n

#endif /* PRECICE_M2N_DISTRIBUTED_COMMUNICATION_HPP_ */
