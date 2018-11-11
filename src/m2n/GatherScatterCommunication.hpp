#pragma once

#include "DistributedCommunication.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"

namespace precice
{
namespace m2n
{

/**
 * @brief Implements DistributedCommunication by using a gathering/scattering methodology.
 * Arrays of data are always gathered and scattered at the master. No direct communication
 * between slaves is used.
 * For more details see m2n/DistributedCommunication.hpp
 */
class GatherScatterCommunication : public DistributedCommunication
{
public:
  GatherScatterCommunication(
      com::PtrCommunication com,
      mesh::PtrMesh         mesh);

  virtual ~GatherScatterCommunication();

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
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  virtual void acceptConnection(
      const std::string &acceptorName,
      const std::string &requesterName);

  /**
   * @brief Requests connection from participant, which has to call acceptConnection().
   *
   * If several connections are going in to a server, the clients have to call this
   * method, while the server has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] nameReuester Name of calling participant.
   */
  virtual void requestConnection(
      const std::string &acceptorName,
      const std::string &requesterName);

  /** same as acceptconnection, but this one does not need vertex distribution
      and instead gets communication map directly from mesh. 
   
   *  This one is used only to create initial communication Map.    
   */
  virtual void acceptPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName);
  
  /** same as requestConnection, but this one does not need vertex distribution
      and instead gets communication map directly from mesh. 
   
   *  This one is used only to create initial communication Map.    
   */
  virtual void requestPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName);


  /** This function should be called by connection accepter to update the vertex list in the 
      communication map, which has be filled previously with demo 
      rank (-1)
  */
  virtual void updateAcceptorCommunicationMap();

  /** This function should be called by connection requester to update the vertex list in the 
      communication map, which has be filled previously with demo 
      rank (-1)
  */
  
  virtual void updateRequesterCommunicationMap();

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  virtual void closeConnection();

  /// Sends an array of double values from all slaves (different for each slave).
  virtual void send(
      double *itemsToSend,
      size_t  size,
      int     valueDimension);

  /// All slaves receive an array of doubles (different for each slave).
  virtual void receive(
      double *itemsToReceive,
      size_t  size,
      int     valueDimension);

  virtual void sendMesh(
    mesh::Mesh &mesh);
  
  /**
   * All ranks receive mesh partition from remote local ranks.
   */
  virtual void receiveMesh(
    mesh::Mesh &mesh);

  /**
   * All ranks Send their local communication maps to connected ranks
   */
  virtual void sendCommunicationMap(
    mesh::Mesh::FeedbackMap &localCommunicationMap);

  /**
   * Each rank revives local communication maps from connected ranks
   */
  virtual void receiveCommunicationMap(
    mesh::Mesh::FeedbackMap &localCommunicationMap) ;

private:
  logging::Logger _log{"m2n::GatherScatterCommunication"};

  /// master to master basic communication
  com::PtrCommunication _com;

  /// Global communication is set up or not
  bool _isConnected;
};

} // namespace m2n
} // namespace precice
