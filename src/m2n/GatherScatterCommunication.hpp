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
      const std::string &requesterName);
  /** 
   *  This method has not been implemented yet.    
   *  @todo: Ideally this should not be here
   */
  virtual void acceptPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName);

  /** 
   *  This method has not been implemented yet.    
   *  @todo: Ideally this should not be here
   */
  virtual void requestPreConnection(
    std::string const &acceptorName,
    std::string const &requesterName);

  /** This function must be called by both accepter and requester to update the vertex list in the 
   *  mapping
  */
  virtual void updateVertexList();
  
  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  void closeConnection() override;

  /// Sends an array of double values from all slaves (different for each slave).
  void send(
      double const *itemsToSend,
      size_t  size,
      int     valueDimension) override;

  /// All slaves receive an array of doubles (different for each slave).
  void receive(
      double *itemsToReceive,
      size_t  size,
      int     valueDimension) override;

   /**
   * @brief Broadcasts an int to connected ranks       
   *        This method has not been implemented yet.    
   *
   * @todo: Ideally this should not be here
   */
  void broadcastSend(const int &itemToSend) override;

  /**
   * @brief Receives an int from a connected rank
   *        This method has not been implemented yet.    
   *
   * @todo: Ideally this should not be here
   */
  void broadcastReceive(std::vector<int> &itemToReceive) override;

  /**
   * @brief All ranks send their mesh partition to remote local  connected ranks.
   *
   * @todo: Ideally this should not be here
   */
  void broadcastSendMesh() override;
  
  /**
   * @brief All ranks receive mesh partition from remote local ranks.
   *
   * @todo: Ideally this should not be here
   */
  void broadcastReceiveMesh() override;

  /**
   *  All ranks Send their local communication maps to connected ranks
   */
  void broadcastSendLCM(
    CommunicationMap &localCommunicationMap) override;

  /*
   *  Each rank revives local communication maps from connected ranks
   */
  void broadcastReceiveLCM(
    CommunicationMap &localCommunicationMap) override;

private:
  logging::Logger _log{"m2n::GatherScatterCommunication"};

  /// master to master basic communication
  com::PtrCommunication _com;

  /// Global communication is set up or not
  bool _isConnected;
};

} // namespace m2n
} // namespace precice
