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

/**TODO
 * @brief Interface for all interprocess communication classes.
 *
 * By default, communication is done within the local communication space. In
 * order to connect to a different communication space, i.e. coupling participant,
 * the methods acceptConnection() and requestConnection() have to called by the
 * two participants which intend to establish a connection. All following
 * communication and process ranking refers to the remote communication space
 * afterwards.
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
   * @brief Sends a std::string to process with given rank.
   */
  virtual void sendMaster (
    const std::string& itemToSend,
    int                rankReceiver ) =0;

  /**
   * @brief Sends an array of integer values.
   */
  virtual void sendMaster (
    int* itemsToSend,
    int  size,
    int  rankReceiver ) =0;

  /**
   * @brief Sends an array of double values.
   */
  virtual void sendMaster (
    double* itemsToSend,
    int     size,
    int     rankReceiver ) =0;

  /**
   * @brief Sends a double to process with given rank.
   */
  virtual void sendMaster (
    double itemToSend,
    int    rankReceiver ) =0;

  /**
   * @brief Sends an int to process with given rank.
   */
  virtual void sendMaster (
    int itemToSend,
    int rankReceiver ) =0;

  /**
   * @brief Sends a bool to process with given rank.
   */
  virtual void sendMaster (
    bool itemToSend,
    int  rankReceiver ) =0;

  /**
   * @brief Receives a std::string from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    std::string& itemToReceive,
    int          rankSender ) =0;

  /**
   * @brief Receives an array of integer values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    int* itemsToReceive,
    int  size,
    int  rankSender ) =0;

  /**
   * @brief Receives an array of double values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    double* itemsToReceive,
    int     size,
    int     rankSender ) =0;

  /**
   * @brief Receives a double from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    double& itemToReceive,
    int     rankSender ) =0;

  /**
   * @brief Receives an int from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    int& itemToReceive,
    int  rankSender ) =0;

  /**
   * @brief Receives a bool from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    bool& itemToReceive,
    int   rankSender ) =0;

  /**
   * @brief Sends an array of double values.
   */
  virtual void sendAll (
    utils::DynVector*   itemsToSend,
    int           size,
    int           rankReceiver,
    mesh::PtrMesh mesh,
    int           valueDimension ) =0;

  virtual void sendAll (
    bool   itemToSend,
    int    rankReceiver) =0;

  virtual void sendAll (
    double itemToSend,
    int    rankReceiver) =0;

  /**
   * @brief Receives an array of double values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual void receiveAll (
    utils::DynVector*   itemsToReceive,
    int           size,
    int           rankSender,
    mesh::PtrMesh mesh,
    int           valueDimension ) =0;

  virtual void receiveAll (
    bool&  itemToReceive,
    int    rankSender ) =0;

  virtual void receiveAll (
    double&  itemToReceive,
    int      rankSender ) =0;



};



}} // namespace precice, m2n

#endif /* PRECICE_M2N_GLOBAL_COMMUNICATION_HPP_ */
