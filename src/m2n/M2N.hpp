// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_M2N_M2N_HPP_
#define PRECICE_M2N_M2N_HPP_

#include "DistributedComFactory.hpp"

#include "com/Communication.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"

#include <map>

namespace precice {
namespace m2n {

/**
 * @brief M2N communication classes for a collection of meshes.
 *
 * TODO
 *
 */
class M2N
{
public:
  using SharedPointer = std::shared_ptr<M2N>;

public:

  M2N( com::Communication::SharedPointer masterCom, DistributedComFactory::SharedPointer distrFactory);

  /**
   * @brief Destructor, empty.
   */
  ~M2N();

  /**
   * @brief Returns true, if a connection to a remote participant has been setup.
   */
  bool isConnected();

  /**
   * @brief Connects to another participant, which has to call requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  void acceptMasterConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester);

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  void requestMasterConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester);


  /**
   * @brief Connects to another participant, which has to call requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  void acceptSlavesConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester);

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  void requestSlavesConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester);

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  void closeConnection();

  /**
   * @brief Get the basic communication between the 2 masters.
   */
  com::Communication::SharedPointer getMasterCommunication();


  void createDistributedCommunication(mesh::PtrMesh mesh);

  void startSendPackage ( int rankReceiver );

  void finishSendPackage();

  /**
   * @brief Starts to receive messages from rankSender.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  int startReceivePackage ( int rankSender );

  void finishReceivePackage();


  /**
   * @brief Sends an array of double values from all slaves (different for each slave).
   */
  void send (
    double* itemsToSend,
    int     size,
    int     meshID,
    int     valueDimension );

  /**
   * @brief The master sends a bool to the other master, for performance reasons, we
   * neglect the gathering and checking step.
   */
  void send (
    bool   itemToSend);

  /**
   * @brief The master sends a double to the other master, for performance reasons, we
   * neglect the gathering and checking step.
   */
  void send (
    double itemToSend);

  /**
   * @brief All slaves receive an array of doubles (different for each slave).
   */
  void receive (
    double* itemsToReceive,
    int     size,
    int     meshID,
    int     valueDimension );

  /**
   * @brief All slaves receive a bool (the same for each slave).
   */
  void receive (
    bool&  itemToReceive);

  /**
   * @brief All slaves receive a double (the same for each slave).
   */
  void receive (
    double&  itemToReceive);

private:

  static tarch::logging::Log _log;

  std::map<int, DistributedCommunication::SharedPointer> _distComs;

  com::Communication::SharedPointer _masterCom;

  DistributedComFactory::SharedPointer _distrFactory;

  bool _isMasterConnected;

  bool _areSlavesConnected;


};



}} // namespace precice, m2n

#endif /* PRECICE_M2N_M2N_HPP_ */
