#pragma once

#include <set>
#include <stddef.h>
#include <string>
#include <vector>

#include "Request.hpp"
#include "boost/range/irange.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "precice/types.hpp"
#include "utils/span.hpp"

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
 *
 * @attention All receive methods, that accept a raw array, expect it to be
 * sized appropriatly. Asynchronous receive methods also expect the vector
 * be sized correctly.
 */
class Communication {

public:
  Communication &operator=(Communication &&) = delete;

  /// Destructor, empty.
  virtual ~Communication()
  {
  }

  /// @name Connection Setup
  /// @{

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
   * @brief Returns a range over all valid remote ranks.
   *
   * @pre A connection to the remote participant has been set up.
   *
   * @see getRemoteCommunicatorSize()
   */
  auto remoteCommunicatorRanks()
  {
    return boost::irange<Rank>(0, static_cast<Rank>(getRemoteCommunicatorSize()));
  }

  /**
   * @brief Accepts connection from another communicator, which has to call requestConnection().
   *
   * Establishes a 1-to-N communication, whereas the acceptor's side is the "1". Contrary to
   * acceptConnectionAsServer(), the other side needs to be a proper communicator with ranks
   * from 0 to N-1. It is not necessary to know this "N" a-priori on the acceptor's side.
   * This communication is used for the 1:1 communication between two master ranks
   * and for the master-slave communication. For the last case, setRankOffset() has to be set.
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   * @param[in] tag Tag for establishing this connection
   * @param[in] acceptorRank Rank of the accpeting process, usually the calling one.
   */
  virtual void acceptConnection(std::string const &acceptorName,
                                std::string const &requesterName,
                                std::string const &tag,
                                int                acceptorRank,
                                int                rankOffset = 0) = 0;

  /**
   * @brief Accepts connection from another communicator, which has to call requestConnectionAsClient().
   *
   * Establishes a 1-to-N communication, whereas the acceptor's side is the "1". Contrary to
   * acceptConnection(), the other side can have arbitrary ranks. However, we need to know its
   * size "N" a-priori.
   * This communication is only used in PointToPointCommunication, i.e. for the M-to-N communication
   * between two participants.
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   * @param[in] tag Tag for establishing this connection
   * @param[in] acceptorRank Rank of accepting server, usually the rank of the current process.
   * @param[in] requesterCommunicatorSize Size of the requester (N)
   */
  virtual void acceptConnectionAsServer(std::string const &acceptorName,
                                        std::string const &requesterName,
                                        std::string const &tag,
                                        int                acceptorRank,
                                        int                requesterCommunicatorSize) = 0;

  /**
   * @brief Connects to another communicator, which has to call acceptConnection().
   *
   * Establishes a 1-to-N communication, whereas the requestor's side is the "N". Contrary to
   * requestConnectionAsClient(), this side needs to be a proper communicator with ranks
   * from 0 to N-1. All ranks need to call this function.
   * This communication is used for the 1:1 communication between two master ranks,
   * and for the master-slave communication.
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   * @param[in] tag Tag for establishing this connection
   * @param[in] requesterRank Rank of the requester (has to go from 0 to N-1)
   * @param[in] requesterCommunicatorSize Size of the requester (N)
   */
  virtual void requestConnection(std::string const &acceptorName,
                                 std::string const &requesterName,
                                 std::string const &tag,
                                 int                requesterRank,
                                 int                requesterCommunicatorSize) = 0;

  /**
   * @brief Connects to another communicator, which has to call acceptConnectionAsServer().
   *
   * Establishes a 1-to-N communication, whereas the requestor's side is the "N". Contrary to
   * requestConnection(), this side can have arbitrary ranks (e.g. 2,3,7). All ranks need to
   * call this function. This communication is only used in PointToPointCommunication, i.e.
   * for the M-to-N communication between two participants.
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to
   * @param[in] tag Tag for establishing this connection
   * @param[in] acceptorRanks Set of ranks that accept a connection
   * @param[in] requesterRank Rank that requests the connection, usually the caller's rank
   */
  virtual void requestConnectionAsClient(std::string const &  acceptorName,
                                         std::string const &  requesterName,
                                         std::string const &  tag,
                                         std::set<int> const &acceptorRanks,
                                         int                  requesterRank) = 0;

  /** Establishes the Master-Slave connection.
   *
   * @param[in] participantName Name of the calling participant.
   * @param[in] tag Tag for establishing this connection
   * @param[in] rank The current rank in the participant 
   * @param[in] size Total size of the participant
   *
   */
  void connectMasterSlaves(std::string const &participantName,
                           std::string const &tag,
                           int                rank,
                           int                size);

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection() = 0;

  /**
   * @brief Prepare environment used to establish the communication.
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  virtual void prepareEstablishment(std::string const &acceptorName,
                                    std::string const &requesterName) {}

  /**
   * @brief Clean-up environment used to establish the communication.
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  virtual void cleanupEstablishment(std::string const &acceptorName,
                                    std::string const &requesterName) {}

  /// @}

  /// @name Reduction
  /// @{

  /// Performs a reduce summation on the rank given by rankMaster
  virtual void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank rankMaster);
  /// Performs a reduce summation on the master, every other rank has to call reduceSum
  virtual void reduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive);

  virtual void reduceSum(int itemToSend, int &itemToReceive, Rank rankMaster);
  virtual void reduceSum(int itemsToSend, int &itemsToReceive);

  virtual void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive, Rank rankMaster);
  virtual void allreduceSum(precice::span<double const> itemsToSend, precice::span<double> itemsToReceive);

  virtual void allreduceSum(double itemToSend, double &itemToReceive, Rank rankMaster);
  virtual void allreduceSum(double itemToSend, double &itemToReceive);

  virtual void allreduceSum(int itemToSend, int &itemToReceive, Rank rankMaster);
  virtual void allreduceSum(int itemToSend, int &itemToReceive);

  /// @}

  /// @name Broadcast
  /// @{

  virtual void broadcast(precice::span<const int> itemsToSend);
  virtual void broadcast(precice::span<int> itemsToReceive, Rank rankBroadcaster);

  virtual void broadcast(int itemToSend);
  virtual void broadcast(int &itemToReceive, Rank rankBroadcaster);

  virtual void broadcast(precice::span<const double> itemsToSend);
  virtual void broadcast(precice::span<double> itemsToReceive, Rank rankBroadcaster);

  virtual void broadcast(double itemToSend);
  virtual void broadcast(double &itemToReceive, Rank rankBroadcaster);

  virtual void broadcast(bool itemToSend);
  virtual void broadcast(bool &itemToReceive, Rank rankBroadcaster);

  virtual void broadcast(std::vector<int> const &v);
  virtual void broadcast(std::vector<int> &v, Rank rankBroadcaster);

  virtual void broadcast(std::vector<double> const &v);
  virtual void broadcast(std::vector<double> &v, Rank rankBroadcaster);

  /// @}

  /// @name Send
  /// @{

  /// Sends a std::string to process with given rank.
  virtual void send(std::string const &itemToSend, Rank rankReceiver) = 0;

  /// Sends an array of integer values.
  virtual void send(precice::span<const int> itemsToSend, Rank rankReceiver) = 0;

  /// Asynchronously sends an array of integer values.
  /// @attention The caller must guarantee that the lifetime of the item extends to the completion of the request!
  virtual PtrRequest aSend(precice::span<const int> itemsToSend, Rank rankReceiver) = 0;

  /// Sends an array of double values.
  virtual void send(precice::span<const double> itemsToSend, Rank rankReceiver) = 0;

  /// Asynchronously sends an array of double values.
  /// @attention The caller must guarantee that the lifetime of the item extends to the completion of the request!
  virtual PtrRequest aSend(precice::span<const double> itemsToSend, Rank rankReceiver) = 0;

  /// @attention The caller must guarantee that the lifetime of the item extends to the completion of the request!
  virtual PtrRequest aSend(std::vector<double> const &itemsToSend, Rank rankReceiver) = 0;

  /// Sends a double to process with given rank.
  virtual void send(double itemToSend, Rank rankReceiver) = 0;

  /// Asynchronously sends a double to process with given rank.
  /// @attention The caller must guarantee that the lifetime of the item extends to the completion of the request!
  virtual PtrRequest aSend(const double &itemToSend, Rank rankReceiver) = 0;

  /// Sends an int to process with given rank.
  virtual void send(int itemToSend, Rank rankReceiver) = 0;

  /// Asynchronously sends an int to process with given rank.
  /// @attention The caller must guarantee that the lifetime of the item extends to the completion of the request!
  virtual PtrRequest aSend(const int &itemToSend, Rank rankReceiver) = 0;

  /// Sends a bool to process with given rank.
  virtual void send(bool itemToSend, Rank rankReceiver) = 0;

  /// Asynchronously sends a bool to process with given rank.
  /// @attention The caller must guarantee that the lifetime of the item extends to the completion of the request!
  virtual PtrRequest aSend(const bool &itemToSend, Rank rankReceiver) = 0;

  /// @}

  /// @name Receive
  /// @{

  /// Receives a std::string from process with given rank.
  virtual void receive(std::string &itemToReceive, Rank rankSender) = 0;

  /// Receives an array of integer values.
  virtual void receive(precice::span<int> itemsToReceive, Rank rankSender) = 0;

  /// Receives an array of double values.
  virtual void receive(precice::span<double> itemsToReceive, Rank rankSender) = 0;

  /// Asynchronously receives an array of double values.
  virtual PtrRequest aReceive(precice::span<double> itemsToReceive, int rankSender) = 0;

  /// Asynchronously receives a vector of double values.
  /*
   * @attention All asynchronous receives methods require the vector to be appropriately sized
   */
  virtual PtrRequest aReceive(std::vector<double> &itemsToReceive, Rank rankSender) = 0;

  /// Receives a double from process with given rank.
  virtual void receive(double &itemToReceive, Rank rankSender) = 0;

  /// Asynchronously receives a double from process with given rank.
  virtual PtrRequest aReceive(double &itemToReceive, Rank rankSender) = 0;

  /// Receives an int from process with given rank.
  virtual void receive(int &itemToReceive, Rank rankSender) = 0;

  /// Asynchronously receives an int from process with given rank.
  virtual PtrRequest aReceive(int &itemToReceive, Rank rankSender) = 0;

  /// Receives a bool from process with given rank.
  virtual void receive(bool &itemToReceive, Rank rankSender) = 0;

  /// Asynchronously receives a bool from process with given rank.
  virtual PtrRequest aReceive(bool &itemToReceive, Rank rankSender) = 0;

  virtual void send(std::vector<int> const &v, Rank rankReceiver) = 0;
  /// Receives an std::vector of ints. The vector will be resized accordingly.
  virtual void receive(std::vector<int> &v, Rank rankSender) = 0;

  virtual void send(std::vector<double> const &v, Rank rankReceiver) = 0;
  /// Receives an std::vector of doubles. The vector will be resized accordingly.
  virtual void receive(std::vector<double> &v, Rank rankSender) = 0;

  /// @}

  /// Set rank offset.
  void setRankOffset(Rank rankOffset)
  {
    _rankOffset = rankOffset;
  }

protected:
  /// Rank offset for masters-slave communication, since ranks are from 0 to size-2
  int _rankOffset = 0;

  bool _isConnected = false;

  /// Adjusts the given rank bases on the _rankOffset
  virtual int adjustRank(Rank rank) const;

private:
  logging::Logger _log{"com::Communication"};
};
} // namespace com
} // namespace precice
