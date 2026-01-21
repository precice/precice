#pragma once

#include <stddef.h>
#include <string>
#include <vector>
#include "DistributedCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice::m2n {

/**
 * @brief Implements DistributedCommunication by using a gathering/scattering methodology.
 * Arrays of data are always gathered and scattered at the primary. No direct communication
 * between secondary ranks is used.
 * For more details see m2n/DistributedCommunication.hpp
 */
class GatherScatterCommunication : public DistributedCommunication {
public:
  GatherScatterCommunication(
      com::PtrCommunication com,
      mesh::PtrMesh         mesh);

  ~GatherScatterCommunication() override;

  /**
   * @brief Returns true, if a connection to a remote participant has been setup.
   */
  bool isConnected() const override;

  /**
   * @brief Accepts connection from participant, which has to call requestConnection().
   *
   * If several connections are going in to a server, the server has to call this
   * method, while the clients have to call requestConnection().
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  void acceptConnection(
      const std::string &acceptorName,
      const std::string &requesterName) override;

  /**
   * @brief Requests connection from participant, which has to call acceptConnection().
   *
   * If several connections are going in to a server, the clients have to call this
   * method, while the server has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] nameReuester Name of calling participant.
   */
  void requestConnection(
      const std::string &acceptorName,
      const std::string &requesterName) override;
  /**
   *  This method has not been implemented yet.
   *  @todo: Ideally this should not be here
   */
  void acceptPreConnection(
      std::string const &acceptorName,
      std::string const &requesterName) override;

  /**
   * @brief Set up a server, so that another participant, can connect using requestConnection().
   *        Call finishAcceptConnection() afterward to respond to requestConnection() calls.
    *
     * @param[in] acceptorName Name of calling participant.
     * @param[in] requesterName Name of remote participant to connect to.
     * @return connectionInfo Information on how to connect to the server that has been set up in this function
     */
    std::string prepareAcceptPreConnection(
        const std::string &acceptorName,
        const std::string &requesterName) override;

  /**
   * @brief Connects to another participant, which has to call requestConnection().
   *        Must not be called before prepareAcceptConnection
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  void finishAcceptPreConnection(
      const std::string &acceptorName,
      const std::string &requesterName) override;

  /**
   *  This method has not been implemented yet.
   *  @todo: Ideally this should not be here
   */
  void requestPreConnection(
      std::string const &acceptorName,
      std::string const &requesterName) override;

  /**
   * @brief Connects to another participant, which has to call acceptPreConnection().
   *        Exchanged vertex list is not included, only connection between ranks
   *        is established.
   *        This version of the function takes a connectionInfo string, that tells
   *        it how to connect to the accepting side.
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   * @param[in] connectionInfo Connection info string containing information on how to connect to accepting side
   */
  void requestPreConnection(
      std::string const                                     &acceptorName,
      std::string const                                     &requesterName,
      com::serialize::SerializedConnectionInfoMap::ConnectionInfoMap const &connectionInfo) override;

  /// Completes the secondary connections for both acceptor and requester by updating the vertex list in _mappings
  void completeSecondaryRanksConnection() override;

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  void closeConnection() override;

  /// Sends an array of double values from all ranks (different for each rank).
  void send(precice::span<double const> itemsToSend, int valueDimension) override;

  /// All ranks receive an array of doubles (different for each rank).
  void receive(precice::span<double> itemsToReceive, int valueDimension) override;

  /// Broadcasts an int to connected ranks on remote participant. Not available for GatherScatterCommunication.
  void broadcastSend(int itemToSend) override;

  /**
   * @brief Receives an int per connected rank on remote participant. Not available for GatherScatterCommunication.
   * @para[out] itemToReceive received ints from remote ranks are stored with the sender rank order
   */
  void broadcastReceiveAll(std::vector<int> &itemToReceive) override;

  /// Broadcasts a mesh to connected ranks on remote participant. Not available for GatherScatterCommunication.
  void broadcastSendMesh() override;

  /// Receive mesh partitions per connected rank on remote participant. Not available for GatherScatterCommunication.
  void broadcastReceiveAllMesh() override;

  /// Scatters a communication map over connected ranks on remote participant. Not available for GatherScatterCommunication.
  void scatterAllCommunicationMap(CommunicationMap &localCommunicationMap) override;

  /// Gathers a communication maps from connected ranks on remote participant. Not available for GatherScatterCommunication.
  void gatherAllCommunicationMap(CommunicationMap &localCommunicationMap) override;

private:
  logging::Logger _log{"m2n::GatherScatterCommunication"};

  /// primary to primary basic communication
  com::PtrCommunication _com;

  /// Global communication is set up or not
  bool _isConnected;
};

} // namespace precice::m2n
