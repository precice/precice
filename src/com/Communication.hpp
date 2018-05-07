#pragma once

#include "Request.hpp"
#include "logging/Logger.hpp"

namespace precice
{
namespace com
{
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
class Communication
{

public:
  /// Destructor, empty.
  virtual ~Communication()
  {
  }

  /// Returns true, if a connection to a remote participant has been setup.
  virtual bool isConnected()
  {
    return _isConnected;
  }

  /**
   * @brief Returns the number of processes in the remote communicator.
   *
   * @pre A connection to the remote participant has been set up.
   */
  virtual size_t getRemoteCommunicatorSize() = 0;

  /**
   * @brief Accepts connection from another communicator, which has to call requestConnection().
   *
   * Establishes a 1-to-N communication, whereas the acceptor's side is the "1". Contrary to
   * acceptConnectionAsServer(), the other side needs to be a proper communicator with ranks
   * from 0 to N-1. It is not necessary to know this "N" a-priori on the acceptor's side.
   * This communication is used for the 1:1 communication between two master ranks, for the
   * 1:N communication to a preCICE server, and for the master-slave communication. For the
   * last case, setRankOffset() has to be set.
   *
   * @param[in] nameAcceptor Name of calling participant.
   * @param[in] nameRequester Name of remote participant to connect to.
   * @param[in] acceptorProcessRank Rank of the acceptor (typically 0)
   * @param[in] acceptorCommunicatorSize Size of the acceptor (has to be 1 assumably)
   */
  virtual void acceptConnection(std::string const &nameAcceptor,
                                std::string const &nameRequester,
                                int                acceptorProcessRank,
                                int                acceptorCommunicatorSize) = 0;

  /**
   * @brief Accepts connection from another communicator, which has to call requestConnectionAsClient().
   *
   * Establishes a 1-to-N communication, whereas the acceptor's side is the "1". Contrary to
   * acceptConnection(), the other side can have arbitrary ranks. However, we need to know its
   * size "N" a-priori.
   * This communication is only used in PointToPointCommunication, i.e. for the M-to-N communication
   * between two participants.
   *
   * @param[in] nameAcceptor Name of calling participant.
   * @param[in] nameRequester Name of remote participant to connect to.
   * @param[in] requesterCommunicatorSize Size of the requestor (N)
   */
  virtual void acceptConnectionAsServer(std::string const &nameAcceptor,
                                        std::string const &nameRequester,
                                        int                requesterCommunicatorSize) = 0;

  /**
   * @brief Connects to another communicator, which has to call acceptConnection().
   *
   * Establishes a 1-to-N communication, whereas the requestor's side is the "N". Contrary to
   * requestConnectionAsClient(), this side needs to be a proper communicator with ranks
   * from 0 to N-1. All ranks need to call this function.
   * This communication is used for the 1:1 communication between two master ranks, for the
   * 1:N communication to a preCICE server, and for the master-slave communication.
   *
   * @param[in] nameAcceptor Name of remote participant to connect to.
   * @param[in] nameRequester Name of calling participant.
   * @param[in] requesterProcessRank Rank of the requestor (has to go from 0 to N-1)
   * @param[in] requesterCommunicatorSize Size of the requestor (N)
   */
  virtual void requestConnection(std::string const &nameAcceptor,
                                 std::string const &nameRequester,
                                 int                requesterProcessRank,
                                 int                requesterCommunicatorSize) = 0;

  /**
   * @brief Connects to another communicator, which has to call acceptConnectionAsServer().
   *
   * Establishes a 1-to-N communication, whereas the requestor's side is the "N". Contrary to
   * requestConnection(), this side can have arbitrary ranks (e.g. 2,3,7). All ranks need to
   * call this function. This communication is only used in PointToPointCommunication, i.e.
   * for the M-to-N communication between two participants.
   *
   * @param[in] nameAcceptor Name of calling participant.
   * @param[in] nameRequester Name of remote participant to connect to.
   */
  virtual int requestConnectionAsClient(std::string const &nameAcceptor,
                                        std::string const &nameRequester) = 0;

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection() = 0;

  /// Performs a reduce summation on the rank given by rankMaster
  virtual void reduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster);

  /// Performs a reduce summation on the master, every other rank has to call reduceSum
  virtual void reduceSum(double *itemsToSend, double *itemsToReceive, int size);

  virtual void reduceSum(int itemToSend, int &itemToReceive, int rankMaster);

  virtual void reduceSum(int itemsToSend, int &itemsToReceive);

  virtual void allreduceSum();

  virtual void allreduceSum(double *itemsToSend, double *itemsToReceive, int size, int rankMaster);

  virtual void allreduceSum(double *itemsToSend, double *itemsToReceive, int size);

  virtual void allreduceSum(double itemToSend, double &itemToReceive, int rankMaster);

  virtual void allreduceSum(double itemToSend, double &itemToReceive);

  virtual void allreduceSum(int itemToSend, int &itemToReceive, int rankMaster);

  virtual void allreduceSum(int itemToSend, int &itemToReceive);

  virtual void broadcast();
  
  virtual void broadcast(int *itemsToSend, int size);

  virtual void broadcast(int *itemsToReceive, int size, int rankBroadcaster);

  virtual void broadcast(int itemToSend);

  virtual void broadcast(int &itemToReceive, int rankBroadcaster);

  virtual void broadcast(double *itemsToSend, int size);

  virtual void broadcast(double *itemsToReceive, int size, int rankBroadcaster);

  virtual void broadcast(double itemToSend);

  virtual void broadcast(double &itemToReceive, int rankBroadcaster);

  virtual void broadcast(bool itemToSend);

  virtual void broadcast(bool &itemToReceive, int rankBroadcaster);

  /// Sends a std::string to process with given rank.
  virtual void send(std::string const &itemToSend, int rankReceiver) = 0;

  /// Sends an array of integer values.
  virtual void send(int *itemsToSend, int size, int rankReceiver) = 0;

  /// Asynchronously sends an array of integer values.
  virtual PtrRequest aSend(int *itemsToSend,
                           int  size,
                           int  rankReceiver) = 0;

  /// Sends an array of double values.
  virtual void send(double *itemsToSend, int size, int rankReceiver) = 0;

  /// Asynchronously sends an array of double values.
  virtual PtrRequest aSend(double *itemsToSend,
                           int     size,
                           int     rankReceiver) = 0;

  /// Sends a double to process with given rank.
  virtual void send(double itemToSend, int rankReceiver) = 0;

  /// Asynchronously sends a double to process with given rank.
  virtual PtrRequest aSend(double *itemToSend,
                           int     rankReceiver) = 0;

  /// Sends an int to process with given rank.
  virtual void send(int itemToSend, int rankReceiver) = 0;

  /// Asynchronously sends an int to process with given rank.
  virtual PtrRequest aSend(int *itemToSend, int rankReceiver) = 0;

  /// Sends a bool to process with given rank.
  virtual void send(bool itemToSend, int rankReceiver) = 0;

  /// Asynchronously sends a bool to process with given rank.
  virtual PtrRequest aSend(bool *itemToSend, int rankReceiver) = 0;

  /// Receives a std::string from process with given rank.
  virtual void receive(std::string &itemToReceive, int rankSender) = 0;

  /// Receives an array of integer values.
  virtual void receive(int *itemsToReceive, int size, int rankSender) = 0;

  /// Asynchronously receives an array of integer values.
  virtual PtrRequest aReceive(int *itemsToReceive,
                              int  size,
                              int  rankSender) = 0;

  /// Receives an array of double values.
  virtual void receive(double *itemsToReceive, int size, int rankSender) = 0;

  /// Asynchronously receives an array of double values.
  virtual PtrRequest aReceive(double *itemsToReceive,
                              int     size,
                              int     rankSender) = 0;

  /// Receives a double from process with given rank.
  virtual void receive(double &itemToReceive, int rankSender) = 0;

  /// Asynchronously receives a double from process with given rank.
  virtual PtrRequest aReceive(double *itemToReceive,
                              int     rankSender) = 0;

  /// Receives an int from process with given rank.
  virtual void receive(int &itemToReceive, int rankSender) = 0;

  /// Asynchronously receives an int from process with given rank.
  virtual PtrRequest aReceive(int *itemToReceive,
                              int  rankSender) = 0;

  /// Receives a bool from process with given rank.
  virtual void receive(bool &itemToReceive, int rankSender) = 0;

  /// Asynchronously receives a bool from process with given rank.
  virtual PtrRequest aReceive(bool *itemToReceive,
                              int   rankSender) = 0;

  /// Set rank offset.
  void setRankOffset(int rankOffset)
  {
    _rankOffset = rankOffset;
  }

  int rank()
  {
    return _rank;
  }

protected:
  int _rank = -1;

  /// Rank offset for masters-slave communication, since ranks are from 0 to size - 2
  int _rankOffset = 0;

  bool _isConnected = false;

private:
  logging::Logger _log{"com::Communication"};
  
};
} // namespace com
} // namespace precice
