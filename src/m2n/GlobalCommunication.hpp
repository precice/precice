// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_M2N_GLOBAL_COMMUNICATION_HPP_
#define PRECICE_M2N_GLOBAL_COMMUNICATION_HPP_

#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "com/Communication.hpp"
#include <string>

namespace precice {
namespace m2n {

/**
 * @brief Interface for all global solver to solver communication classes.
 *
 * By default, communication is done within the local communication space. In
 * order to connect to a different communication space, i.e. coupling participant,
 * the methods acceptConnection() and requestConnection() have to called by the
 * two participants which intend to establish a connection. All following
 * communication and process ranking refers to the remote communication space
 * afterwards.
 *
 * If both solver are run in serial, this interface falls back to a simple wrapper
 * of the normal com/Communication.hpp . If one is run in parallel, the parallel
 * communication is defined here. The basic communication (e.g. Sockets or MPI) is
 * still handled in com/Communication.hpp .
 *
 * The class offers methods to communicate between the 2 master processes, send/receiveMaster
 * and methods to communicate between the slaves, send/ReceiveAll. This can either
 * mean that data is communicated point2point, in case of arrays, or that single values
 * are broadcasted.
 *
 */
class GlobalCommunication
{
public:

  /**
   * @brief Destructor, empty.
   */
  virtual ~GlobalCommunication() {}

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
    const std::string& nameRequester,
    int                acceptorProcessRank,
    int                acceptorCommunicatorSize ) =0;

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  virtual void requestConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester,
    int                requesterProcessRank,
    int                requesterCommunicatorSize ) =0;

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection() =0;

  /**
   * @brief Get the basic communication between the 2 masters.
   */
  virtual com::PtrCommunication getMasterCommunication() =0;

  virtual void startSendPackage ( int rankReceiver ) =0;

  virtual void finishSendPackage() =0;

  /**
   * @brief Starts to receive messages from rankSender.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int startReceivePackage ( int rankSender ) =0;

  virtual void finishReceivePackage() =0;

  /**
   * @brief Sends a std::string from master to master.
   */
  virtual void sendMaster (
    const std::string& itemToSend,
    int                rankReceiver ) =0;

  /**
   * @brief Sends an array of integer values from master to master.
   */
  virtual void sendMaster (
    int* itemsToSend,
    int  size,
    int  rankReceiver ) =0;

  /**
   * @brief Sends an array of double values from master to master.
   */
  virtual void sendMaster (
    double* itemsToSend,
    int     size,
    int     rankReceiver ) =0;

  /**
   * @brief Sends a double from master to master .
   */
  virtual void sendMaster (
    double itemToSend,
    int    rankReceiver ) =0;

  /**
   * @brief Sends an int from master to master.
   */
  virtual void sendMaster (
    int itemToSend,
    int rankReceiver ) =0;

  /**
   * @brief Sends a bool from master to master.
   */
  virtual void sendMaster (
    bool itemToSend,
    int  rankReceiver ) =0;

  /**
   * @brief Receives a std::string from master to master.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    std::string& itemToReceive,
    int          rankSender ) =0;

  /**
   * @brief Receives an array of integer values from master to master.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    int* itemsToReceive,
    int  size,
    int  rankSender ) =0;

  /**
   * @brief Receives an array of double values from master to master.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    double* itemsToReceive,
    int     size,
    int     rankSender ) =0;

  /**
   * @brief Receives a double from process with given rank from master to master.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    double& itemToReceive,
    int     rankSender ) =0;

  /**
   * @brief Receives an int from process with given rank from master to master.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    int& itemToReceive,
    int  rankSender ) =0;

  /**
   * @brief Receives a bool from process with given rank from master to master.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    bool& itemToReceive,
    int   rankSender ) =0;

  /**
   * @brief Sends an array of double values from all slaves (different for each slave).
   */
  virtual void sendAll (
    utils::DynVector*   itemsToSend,
    int           size,
    int           rankReceiver,
    mesh::PtrMesh mesh,
    int           valueDimension ) =0;

  /**
   * @brief The master sends a bool to the other master, for performance reasons, we
   * neglect the gathering and checking step.
   */
  virtual void sendAll (
    bool   itemToSend,
    int    rankReceiver) =0;

  /**
   * @brief The master sends a double to the other master, for performance reasons, we
   * neglect the gathering and checking step.
   */
  virtual void sendAll (
    double itemToSend,
    int    rankReceiver) =0;

  /**
   * @brief All slaves receive an array of doubles (different for each slave).
   */
  virtual void receiveAll (
    utils::DynVector*   itemsToReceive,
    int           size,
    int           rankSender,
    mesh::PtrMesh mesh,
    int           valueDimension ) =0;

  /**
   * @brief All slaves receive a bool (the same for each slave).
   */
  virtual void receiveAll (
    bool&  itemToReceive,
    int    rankSender ) =0;

  /**
   * @brief All slaves receive a double (the same for each slave).
   */
  virtual void receiveAll (
    double&  itemToReceive,
    int      rankSender ) =0;



};



}} // namespace precice, m2n

#endif /* PRECICE_M2N_GLOBAL_COMMUNICATION_HPP_ */
