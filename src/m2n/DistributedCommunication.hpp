#pragma once

#include "mesh/SharedPointer.hpp"
#include <map>
#include <vector>

namespace precice
{
namespace m2n
{

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
class DistributedCommunication
{
public:
  using SharedPointer = std::shared_ptr<DistributedCommunication>;

  explicit DistributedCommunication(mesh::PtrMesh mesh)
      : _mesh(mesh)
  {}

  /// Destructor, empty.
  virtual ~DistributedCommunication() {}

  /// Returns true, if a connection to a remote participant has been setup.
  virtual bool isConnected() = 0;

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

  /** same as acceptconnection, but this one does not need vertex distribution
      and instead gets communication map directly from mesh. 
   
   *  This one is used only to create initial communication Map.    
   */
  virtual void acceptPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName) = 0;
  
  /** same as requestConnection, but this one does not need vertex distribution
      and instead gets communication map directly from mesh. 
   
   *  This one is used only to create initial communication Map.    
   */
  virtual void requestPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName) = 0;

  /** This function should be called by both accepter and requester to update the vertex list in the 
   *  mapping
  */
  virtual void updateVertexList() = 0;
  
  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection() = 0;

  /// Sends an array of double values from all slaves (different for each slave).
  virtual void send(
      double *itemsToSend,
      size_t  size,
      int     valueDimension) = 0;

  /// All slaves receive an array of doubles (different for each slave).
  virtual void receive(
      double *itemsToReceive,
      size_t  size,
      int     valueDimension) = 0;

  /**
   * @brief Sends a double to connected ranks       
   */
  virtual void broadcastSend(double &itemToSend) = 0;

  /**
   * @brief Receives a double from a connected rank
   */
  virtual void broadcastReceive(double &itemToReceive) = 0;

  /**
   * All ranks send their mesh partition to remote local  connected ranks.
   */
  virtual void broadcastSendMesh(
    mesh::Mesh &mesh) = 0;
  
  /**
   * All ranks receive mesh partition from remote local ranks.
   */
  virtual void broadcastReceiveMesh(
    mesh::Mesh &mesh) = 0;

  /**
   *  All ranks Send their local communication maps to connected ranks
   */
  virtual void broadcastSendLCM(
    std::map<int, std::vector<int>> &localCommunicationMap)=0;

  /*
   *  Each rank revives local communication maps from connected ranks
   */
  virtual void broadcastReceiveLCM(
    std::map<int, std::vector<int>> &localCommunicationMap)=0 ;
  
protected:
  /**
   * @brief mesh that dictates the distribution of this mapping
   *
   * @todo maybe change this directly to vertexDistribution
   */
  mesh::PtrMesh _mesh;
};

} // namespace m2n
} // namespace precice
