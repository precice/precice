#pragma once

#include <list>
#include <vector>
#include <string>

#include "DistributedCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Mesh.hpp"

namespace precice
{
namespace m2n
{
/**
 * @brief Point-to-point communication implementation of DistributedCommunication.
 *
 * Direct communication of local data subsets is performed between processes of
 * coupled participants. The two supported implementations of direct
 * communication are SocketCommunication and MPIPortsCommunication which can be
 * supplied via their corresponding instantiation factories
 * SocketCommunicationFactory and MPIPortsCommunicationFactory.
 *
 * For the detailed implementation documentation refer to PointToPointCommunication.cpp.
 */
class PointToPointCommunication : public DistributedCommunication
{
public:
  PointToPointCommunication(com::PtrCommunicationFactory communicationFactory,
                            mesh::PtrMesh                mesh);

  ~PointToPointCommunication() override;

  /// Returns true, if a connection to a remote participant has been established.
  bool isConnected() const override;

  /**
   * @brief Accepts connection from participant, which has to call
   *        requestConnection().
   *
   * @param[in] acceptorName  Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  void acceptConnection(std::string const &acceptorName,
                        std::string const &requesterName) override;

  /**
   * @brief Requests connection from participant, which has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  void requestConnection(std::string const &acceptorName,
                         std::string const &requesterName) override;

  /** same as acceptconnection, but this one does not need vertex distribution
      and instead gets connected ranks directly from mesh. 
   
   *  This one is used only to create initial communication Map.    
   */
  virtual void acceptPreConnection(std::string const &acceptorName,
                                   std::string const &requesterName);
  
  /** same as requestConnection, but this one does not need vertex distribution
      and instead gets connected ranks directly from mesh. 
   
   *  This one is used only to create initial communication Map.    
   */
  virtual void requestPreConnection(std::string const &acceptorName,
                                    std::string const &requesterName);

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  void closeConnection() override;

  /**
   * @brief Sends a subset of local double values corresponding to local indices
   *        deduced from the current and remote vertex distributions.
   */
  void send(double const *itemsToSend, size_t size, int valueDimension = 1) override;

  /**
   * @brief Receives a subset of local double values corresponding to local
   *        indices deduced from the current and remote vertex distributions.
   */
  void receive(double *itemsToReceive,
               size_t  size,
               int     valueDimension = 1) override;

  /**
   * @brief Broadcasts a double to connected ranks       
   */
  virtual void broadcastSend(const double &itemToSend);

  /**
   * @brief Receives a double from a connected rank
   */
  virtual void broadcastReceive(double &itemToReceive);

  /**
   * All ranks send their mesh partition to remote local  connected ranks.
   */
  virtual void broadcastSendMesh(mesh::Mesh &mesh);
  
  /**
   * All ranks receive mesh partition from remote local ranks.
   */
  virtual void broadcastReceiveMesh(mesh::Mesh &mesh);

  /**
   *  All ranks Send their local communication maps to connected ranks
   */
  virtual void broadcastSendLCM(
    std::map<int, std::vector<int>> &localCommunicationMap);

  /**
   *  Each rank revives local communication maps from connected ranks
   */
  virtual void broadcastReceiveLCM(
    std::map<int, std::vector<int>> &localCommunicationMap);

private:
  logging::Logger _log{"m2n::PointToPointCommunication"};

  /// Checks all stored requests for completion and removes associated buffers
  /**
   * @param[in] blocking False means that the function returns, even when there are requests left.
   */  
  void checkBufferedRequests(bool blocking);
  
  com::PtrCommunicationFactory _communicationFactory;

  /// Communication class used for this PointToPointCommunication
  /**
   * A Communication object represents all connections to all ranks made by this P2P instance.
   **/
  com::PtrCommunication _communication;
  
  /**
   * @brief Defines mapping between:
   *        1. global remote process rank;
   *        2. local data indices, which define a subset of local (for process
   *           rank in the current participant) data to be communicated between
   *           the current process rank and the remote process rank;
   *        3. Request holding information about pending communication
   *        4. Appropriately sized buffer to receive elements
   */
  struct Mapping {
    int                 remoteRank;
    std::vector<int>    indices;
    com::PtrRequest     request;
    std::vector<double> recvBuffer;
  };

  /**
   * @brief Local (for process rank in the current participant) vector of
   *        mappings (one to service each point-to-point connection).
   */
  std::vector<Mapping> _mappings;

   /**
   * @brief this data structure is used to store m2n communication information for the 1 step of 
   *        bounding box initialization. It stores:
   *        1. global remote process rank;
   *        2. communication object (provides point-to-point communication routines).
   *        3. Appropriatly sized buffer to receive elements
   */
  struct ConnectionData {
    int                   remoteRank;
    com::PtrCommunication communication;
    com::PtrRequest       request;
  };

  /**
   * @brief Local (for process rank in the current participant) vector of
   *        ConnectionData (one to service each point-to-point connection).
   */
  std::vector<ConnectionData> _connectionDataVector;

  bool _isConnected = false;

  std::list<std::pair<std::shared_ptr<com::Request>,
                      std::shared_ptr<std::vector<double>>>> bufferedRequests;

};
} // namespace m2n
} // namespace precice
