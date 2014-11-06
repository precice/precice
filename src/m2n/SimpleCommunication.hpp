// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_M2N_SIMPLE_COMMUNICATION_HPP_
#define PRECICE_M2N_SIMPLE_COMMUNICATION_HPP_

#include "GlobalCommunication.hpp"
#include "com/Communication.hpp"
#include "tarch/logging/Log.h"
#include "com/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"


namespace precice {
namespace m2n {

/** TODO
 * @brief Implements Communication by using sockets.
 */
class SimpleCommunication : public GlobalCommunication
{
public:

  /**
   * @brief Constructor.
   */
  SimpleCommunication (
     com::PtrCommunication com);

  /**
   * @brief Destructor.
   */
  virtual ~SimpleCommunication();

  /**
   * @brief Returns true, if a connection to a remote participant has been setup.
   */
  virtual bool isConnected();

  /**
   * @brief Accepts connection from participant, which has to call requestConnection().
   *
   * If several connections are going in to a server, the server has to call this
   * method, while the clients have to call requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  virtual void acceptConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester,
    int                acceptorProcessRank,
    int                acceptorCommunicatorSize );

  /**
   * @brief Requests connection from participant, which has to call acceptConnection().
   *
   * If several connections are going in to a server, the clients have to call this
   * method, while the server has to call acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  virtual void requestConnection (
    const std::string& nameAcceptor,
    const std::string& nameRequester,
    int                requesterProcessRank,
    int                requesterCommunicatorSize );

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection();

  virtual com::PtrCommunication getMasterCommunication();

  /**
   * @brief Is empty.
   */
  virtual void startSendPackage ( int rankReceiver );

  /**
   * @brief Is empty.
   */
  virtual void finishSendPackage();

  /**
   * @brief Just returns rank of sender.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int startReceivePackage ( int rankSender );

  /**
   * @brief Is empty.
   */
  virtual void finishReceivePackage();

  /**
   * @brief Sends a std::string to process with given rank.
   */
  virtual void sendMaster (
    const std::string& itemToSend,
    int                rankReceiver );

  /**
   * @brief Sends an array of integer values.
   */
  virtual void sendMaster (
    int* itemsToSend,
    int  size,
    int  rankReceiver );

  /**
   * @brief Sends an array of double values.
   */
  virtual void sendMaster (
    double* itemsToSend,
    int     size,
    int     rankReceiver );

  /**
   * @brief Sends a double to process with given rank.
   */
  virtual void sendMaster (
    double itemToSend,
    int    rankReceiver );

  /**
   * @brief Sends an int to process with given rank.
   */
  virtual void sendMaster (
    int itemToSend,
    int rankReceiver );

  /**
   * @brief Sends a bool to process with given rank.
   */
  virtual void sendMaster (
    bool itemToSend,
    int  rankReceiver );

  /**
   * @brief Receives a std::string from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    std::string& itemToReceive,
    int          rankSender );

  /**
   * @brief Receives an array of integer values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    int* itemsToReceive,
    int  size,
    int  rankSender );

  /**
   * @brief Receives an array of double values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    double* itemsToReceive,
    int     size,
    int     rankSender );

  /**
   * @brief Receives a double from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    double& itemToReceive,
    int     rankSender );

  /**
   * @brief Receives an int from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    int& itemToReceive,
    int  rankSender );

  /**
   * @brief Receives a bool from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receiveMaster (
    bool& itemToReceive,
    int   rankSender );

  /**
   * @brief Sends an array of double values.
   */
  virtual void sendAll (
    double*       itemsToSend,
    int           size,
    int           rankReceiver,
    mesh::PtrMesh mesh,
    int           valueDimension);


  /**
   * @brief Receives an array of double values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual void receiveAll (
    double*       itemsToReceive,
    int           size,
    int           rankSender,
    mesh::PtrMesh mesh,
    int           valueDimension);

private:

  static tarch::logging::Log _log;

  com::PtrCommunication _com;
};

}} // namespace precice, m2n

#endif /* PRECICE_M2N_SIMPLE_COMMUNICATION_HPP_ */
