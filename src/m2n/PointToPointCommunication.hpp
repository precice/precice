#pragma once

#include "DistributedCommunication.hpp"
#include <list>
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

  virtual ~PointToPointCommunication();

  /// Returns true, if a connection to a remote participant has been established.
  virtual bool isConnected();

  /**
   * @brief Accepts connection from participant, which has to call
   *        requestConnection().
   *
   * @param[in] nameAcceptor  Name of calling participant.
   * @param[in] nameRequester Name of remote participant to connect to.
   */
  virtual void acceptConnection(std::string const &nameAcceptor,
                                std::string const &nameRequester);

  /**
   * @brief Requests connection from participant, which has to call acceptConnection().
   *
   * @param[in] nameAcceptor Name of remote participant to connect to.
   * @param[in] nameRequester Name of calling participant.
   */
  virtual void requestConnection(std::string const &nameAcceptor,
                                 std::string const &nameRequester);

  /** same as acceptconnection, but this one does not need vertex distribution
      and instead gets communication map directly from mesh. 
   
   *  This one is used only to create initial communication Map.    
   */
  virtual void acceptPreConnection(
    std::string const &nameAcceptor,
    std::string const &nameRequester);
  
  /** same as requestConnection, but this one does not need vertex distribution
      and instead gets communication map directly from mesh. 
   
   *  This one is used only to create initial communication Map.    
   */
  virtual void requestPreConnection(
    std::string const &nameAcceptor,
    std::string const &nameRequester);

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection();

  /**
   * @brief Sends a subset of local double values corresponding to local indices
   *        deduced from the current and remote vertex distributions.
   */
  virtual void send(double *itemsToSend, size_t size, int valueDimension = 1);

  /**
   * @brief Receives a subset of local double values corresponding to local
   *        indices deduced from the current and remote vertex distributions.
   */
  virtual void receive(double *itemsToReceive,
                       size_t  size,
                       int     valueDimension = 1);

  /// All ranks Send their partition to remote local ranks.
  virtual void sendMesh(
    mesh::Mesh &mesh);
  
  /// All ranks receive mesh partition from remote local ranks.
  virtual void receiveMesh(
    mesh::Mesh &mesh);

  /// All ranks Send their local communication maps to connected ranks
  virtual void sendCommunicationMap(
    mesh::Mesh::FeedbackMap &localCommunicationMap);

  /// Each rank revives local communication maps from connected ranks
  virtual void receiveCommunicationMap(
    mesh::Mesh::FeedbackMap &localCommunicationMap) ;

private:
  logging::Logger _log{"m2n::PointToPointCommunication"};

  /// Checks all stored requests for completion and removes associated buffers
  /**
   * @param[in] blocking False means that the function returns, even when there are requests left.
   */  
  void checkBufferedRequests(bool blocking);
  
  com::PtrCommunicationFactory _communicationFactory;

  /**
   * @brief Defines mapping between:
   *        1. local (to the current process) remote process rank;
   *        2. global remote process rank;
   *        3. local data indices, which define a subset of local (for process
   *           rank in the current participant) data to be communicated between
   *           the current process rank and the remote process rank;
   *        4. communication object (provides point-to-point communication routines).
   *        5. Appropriatly sized buffer to receive elements
   */
  struct Mapping {
    int                   localRemoteRank;
    int                   globalRemoteRank;
    std::vector<int>      indices;
    com::PtrCommunication communication;
    com::PtrRequest       request;
    std::vector<double>   recvBuffer;
  };

  /**
   * @brief Local (for process rank in the current participant) vector of
   *        mappings (one to service each point-to-point connection).
   */
  std::vector<Mapping> _mappings;

  bool _isConnected = false;

  std::map<int, std::vector<int>> _localCommunicationMap;

  std::list<std::pair<std::shared_ptr<com::Request>,
                      std::shared_ptr<std::vector<double>>>> bufferedRequests;

};
} // namespace m2n
} // namespace precice
