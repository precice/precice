// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at
// http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_COMMUNICATION_HPP_
#define PRECICE_COM_COMMUNICATION_HPP_

#include "Request.hpp"

#include <memory>
#include "tarch/logging/Log.h"

namespace precice {
namespace com {
/**
 * @brief Interface for all interprocess communication classes.
 *
 * By default, communication is done within the local communication space. In
 * order to connect to a different communication space, i.e. coupling
 * participant, the methods acceptConnection() and requestConnection() have to
 * called by the two participants which intend to establish a connection. All
 * following communication and process ranking refers to the remote
 * communication space afterwards.
 *
 * Sending methods prefixed with `a' are asynchronous. It means that they return
 * immediately, even though either the actual sending might have not been
 * started yet or data from user buffer that is being supplied has not been
 * safely stored away (to system buffer) yet. This implies that user buffer
 * cannot be immediately reused (for writing) after asynchronous sending call
 * returns. However, a special corresponding "request" object is returned by all
 * asynchronous sending methods, which could be further used in future in order
 * to properly wait (block execution) until either the corresponding sending
 * request has truly finished or data from user buffer that is being supplied
 * has been safely stored away (to system buffer) and, thus, can be reused (for
 * writing).
 *
 * The main benefit from asynchronous sending methods is their deterministic
 * behavior, i.e. the guarantee that they can never block the execution (return
 * immediately). The two typical scenarios where this comes handy are:
 * 1. Robustness --- allows one to write deadlock-prone communication algorithms
 *    (see point-to-point communication).
 * 2. Efficiency --- allows one to perform some useful computational work while
 *    sending is happening on the background.
 */
class Communication {
public:
  using SharedPointer = std::shared_ptr<Communication>;

public:
  enum {
    // @brief Can be used instead of specific rankSender number in receive
    // calls.
    ANY_SENDER = -1
  };

  Communication() : _rank(-1), _rankOffset(0) {
  }

  /**
   * @brief Destructor, empty.
   */
  virtual ~Communication() {
  }

  /**
   * @brief Returns true, if a connection to a remote participant has been
   * setup.
   */
  virtual bool isConnected() = 0;

  /**
   * @brief Returns the number of processes in the remote communicator.
   *
   * Precondition: a connection to the remote participant has been setup.
   */
  virtual int getRemoteCommunicatorSize() = 0;

  /**
   * @brief Connects to another participant, which has to call
   * requestConnection().
   *
   * @param nameAcceptor [IN] Name of calling participant.
   * @param nameRequester [IN] Name of remote participant to connect to.
   */
  virtual void acceptConnection(std::string const& nameAcceptor,
                                std::string const& nameRequester,
                                int acceptorProcessRank,
                                int acceptorCommunicatorSize) = 0;

  virtual void acceptConnectionAsServer(std::string const& nameAcceptor,
                                        std::string const& nameRequester,
                                        int requesterCommunicatorSize) = 0;

  /**
   * @brief Connects to another participant, which has to call
   * acceptConnection().
   *
   * @param nameAcceptor [IN] Name of remote participant to connect to.
   * @param nameReuester [IN] Name of calling participant.
   */
  virtual void requestConnection(std::string const& nameAcceptor,
                                 std::string const& nameRequester,
                                 int requesterProcessRank,
                                 int requesterCommunicatorSize) = 0;

  virtual int requestConnectionAsClient(std::string const& nameAcceptor,
                                        std::string const& nameRequester) = 0;

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection() = 0;

  virtual void startSendPackage(int rankReceiver) = 0;

  virtual void finishSendPackage() = 0;

  /**
   * @brief Starts to receive messages from rankSender.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int startReceivePackage(int rankSender) = 0;

  virtual void finishReceivePackage() = 0;

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

  /**
   * @brief Sends a std::string to process with given rank.
   */
  virtual void send(std::string const& itemToSend, int rankReceiver) = 0;

  /**
   * @brief Sends an array of integer values.
   */
  virtual void send(int* itemsToSend, int size, int rankReceiver) = 0;

  /**
   * @brief Asynchronously sends an array of integer values.
   */
  virtual Request::SharedPointer aSend(int* itemsToSend,
                                       int size,
                                       int rankReceiver) = 0;

  /**
   * @brief Sends an array of double values.
   */
  virtual void send(double* itemsToSend, int size, int rankReceiver) = 0;

  /**
   * @brief Asynchronously sends an array of double values.
   */
  virtual Request::SharedPointer aSend(double* itemsToSend,
                                       int size,
                                       int rankReceiver) = 0;

  /**
   * @brief Sends a double to process with given rank.
   */
  virtual void send(double itemToSend, int rankReceiver) = 0;

  /**
   * @brief Asynchronously sends a double to process with given rank.
   */
  virtual Request::SharedPointer aSend(double* itemToSend,
                                       int rankReceiver) = 0;

  /**
   * @brief Sends an int to process with given rank.
   */
  virtual void send(int itemToSend, int rankReceiver) = 0;

  /**
   * @brief Asynchronously sends an int to process with given rank.
   */
  virtual Request::SharedPointer aSend(int* itemToSend, int rankReceiver) = 0;

  /**
   * @brief Sends a bool to process with given rank.
   */
  virtual void send(bool itemToSend, int rankReceiver) = 0;

  /**
   * @brief Asynchronously sends a bool to process with given rank.
   */
  virtual Request::SharedPointer aSend(bool* itemToSend, int rankReceiver) = 0;

  /**
   * @brief Receives a std::string from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(std::string& itemToReceive, int rankSender) = 0;

  /**
   * @brief Receives an array of integer values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(int* itemsToReceive, int size, int rankSender) = 0;

  /**
   * @brief Receives an array of double values.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(double* itemsToReceive, int size, int rankSender) = 0;

  /**
   * @brief Receives a double from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(double& itemToReceive, int rankSender) = 0;

  /**
   * @brief Receives an int from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(int& itemToReceive, int rankSender) = 0;

  /**
   * @brief Receives a bool from process with given rank.
   *
   * @return Rank of sender, which is useful when ANY_SENDER is used.
   */
  virtual int receive(bool& itemToReceive, int rankSender) = 0;

  /**
   * @brief Set rank offset.
   */
  void
  setRankOffset(int rankOffset) {
    _rankOffset = rankOffset;
  }

  int
  rank() {
    return _rank;
  }

protected:
  int _rank;

  /**
   * @brief Rank offset for masters-slave communication, since ranks are from 0
   * to size - 2
   */
  int _rankOffset;

private:
  // @brief Logging device.
  static tarch::logging::Log _log;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_COMMUNICATION_HPP_ */
