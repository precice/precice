#pragma once

#include <stddef.h>
#include <string>
#include <vector>
#include "DistributedCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace m2n {

/**
 * @brief Implements DistributedCommunication by using a gathering/scattering methodology.
 * Arrays of data are always gathered and scattered at the master. No direct communication
 * between slaves is used.
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
   *  This method has not been implemented yet.
   *  @todo: Ideally this should not be here
   */
  void requestPreConnection(
      std::string const &acceptorName,
      std::string const &requesterName) override;

  /// Completes the slaves connections for both acceptor and requester by updating the vertex list in _mappings
  void completeSlavesConnection() override;

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  void closeConnection() override;

  /// Sends an array of double values from all slaves (different for each slave).
  void send(
      double const *itemsToSend,
      size_t        size,
      int           valueDimension) override;

  /// All slaves receive an array of doubles (different for each slave).
  void receive(
      double *itemsToReceive,
      size_t  size,
      int     valueDimension) override;

  /// Broadcasts an int to connected ranks on remote participant. Not available for GatherScatterCommunication.
  void broadcastSend(const int &itemToSend) override;

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

  /// master to master basic communication
  com::PtrCommunication _com;

  /// Global communication is set up or not
  bool _isConnected;
};

} // namespace m2n
} // namespace precice
