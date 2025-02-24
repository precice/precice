#pragma once

#include <map>
#include <vector>
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "precice/span.hpp"

namespace precice::m2n {

/**
 * @brief Interface for all distributed solver to solver communication classes.
 *
 *
 * By default, communication is done within the local communication space. In
 * order to connect to a different communication space, i.e. coupling participant,
 * the methods acceptConnection() and requestConnection() have to called by the
 * two participants which intend to establish a connection. All following
 * communication and process ranking refers to the remote communication space
 * afterwards.
 *
 * This interface organizes the communication between 2 distributed participants.
 * The core communication (e.g. Sockets or MPI) is still handled in com/Communication.hpp .
 * The information on how the mesh is distributed can be accessed through the member variable
 * _mesh.
 *
 * The class offers methods to communicate between all processors. This can either
 * mean that data is communicated in a distributed way, in case of arrays, or that single values
 * are broadcasted.
 *
 */
class DistributedCommunication {
public:
  using SharedPointer = std::shared_ptr<DistributedCommunication>;

  explicit DistributedCommunication(mesh::PtrMesh mesh)
      : _mesh(mesh)
  {
  }

  /// Destructor, empty.
  virtual ~DistributedCommunication() = default;

  /// Returns true, if a connection to a remote participant has been setup.
  virtual bool isConnected() const = 0;

  /**
   * @brief Connects to another participant, which has to call requestConnection().
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  virtual void acceptConnection(
      const std::string &acceptorName,
      const std::string &requesterName) = 0;

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  virtual void requestConnection(
      const std::string &acceptorName,
      const std::string &requesterName) = 0;

  /**
   * @brief Connects to another participant, which has to call requestPreConnection().
   *        Exchanged vertex list is not included, only connection between ranks
   *        is established.
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  virtual void acceptPreConnection(
      std::string const &acceptorName,
      std::string const &requesterName) = 0;

  /**
   * @brief Connects to another participant, which has to call acceptPreConnection().
   *        Exchanged vertex list is not included, only connection between ranks
   *        is established.
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  virtual void requestPreConnection(
      std::string const &acceptorName,
      std::string const &requesterName) = 0;

  /* @brief Completes the secondary connections for both acceptor and requester by updating
   * the vertex list in _mappings.
   *
   * @pre acceptPreConnection or requestPreConnection has been called.
   * @pre accept/requestConnection has not been called
   * @post the m2n connection is ready to use
   */
  virtual void completeSecondaryRanksConnection() = 0;

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection() = 0;

  /// Sends an array of double values from all ranks (different for each rank).
  virtual void send(precice::span<double const> itemsToSend, int valueDimension) = 0;

  /// All ranks receive an array of doubles (different for each rank).
  virtual void receive(precice::span<double> itemsToReceive, int valueDimension) = 0;

  /*
   * A mapping from remote local ranks to the IDs that must be communicated
   */
  using CommunicationMap = std::map<int, std::vector<int>>;

  /// Broadcasts an int to connected ranks on remote participant
  virtual void broadcastSend(int itemToSend) = 0;

  /**
   * @brief Receives an int per connected rank on remote participant
   * @para[out] itemToReceive received ints from remote ranks are stored with the sender rank order
   */
  virtual void broadcastReceiveAll(std::vector<int> &itemToReceive) = 0;

  /// Broadcasts a mesh to connected ranks on remote participant
  virtual void broadcastSendMesh() = 0;

  /// Receive mesh partitions per connected rank on remote participant
  virtual void broadcastReceiveAllMesh() = 0;

  /// Scatters a communication map over connected ranks on remote participant
  virtual void scatterAllCommunicationMap(CommunicationMap &localCommunicationMap) = 0;

  /// Gathers a communication maps from connected ranks on remote participant
  virtual void gatherAllCommunicationMap(CommunicationMap &localCommunicationMap) = 0;

protected:
  /**
   * @brief mesh that dictates the distribution of this mapping
   *
   * @todo maybe change this directly to vertexDistribution
   */
  mesh::PtrMesh _mesh;
};

} // namespace precice::m2n
