// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_MPIASCODTCOMMUNICATION_HPP_
#define PRECICE_COM_MPIASCODTCOMMUNICATION_HPP_

#include "mpi.h"


/**
 * @brief Provides implementation for ascodt MPI point-to-point communication.
 */
class MPIAscodtCommunication
{
public:

  /**
   * @brief Constructor
   */
  MPIAscodtCommunication();

  /**
   * @brief Destructor, empty.
   */
  ~MPIAscodtCommunication();

  /**
   * @brief Returns true, if a connection to a remote participant has been setup.
   */
  bool isConnected();

  /**
   * @brief See precice::com::Communication::acceptConnection().
   */
  void acceptConnection (
      int* vertexTable,
      int* adressTable);

  /**
   * @brief See precice::com::Communication::requestConnection().
   */
  void requestConnection (
      int* vertexTable,
      int* adressTable);

  /**
   * @brief See precice::com::Communication::closeConnection().
   */
  void closeConnection();

  /**
   * @brief Sends an array of double values.
   */
  void send (
    double* itemsToSend,
    int     size);

  /**
   * @brief Receives an array of double values.
   */
  void receive (
    double* itemsToReceive,
    int     size);

};


#endif /* PRECICE_COM_MPIASCODTCOMMUNICATION_HPP_ */
