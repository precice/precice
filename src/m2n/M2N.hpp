#pragma once

#include "DistributedComFactory.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"
#include <map>

namespace precice
{
namespace m2n
{

/**
 * @brief M2N communication class.
 * This layer is necessary since communication between two participants can be working via several meshes,
 * each possibly with a different decomposition. In principle, this class is only a map from meshes to DistributedCommunications
 *
 */
class M2N
{
public:
  M2N(com::PtrCommunication masterCom, DistributedComFactory::SharedPointer distrFactory);

  /// Destructor, empty.
  ~M2N();

  /// Returns true, if a connection to a remote participant has been setup.
  bool isConnected();

  /**
   * @brief Connects to another participant, which has to call requestConnection().
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  void acceptMasterConnection(const std::string &acceptorName,
                              const std::string &requesterName);

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  void requestMasterConnection(const std::string &acceptorName,
                               const std::string &requesterName);

  /**
   * @brief Connects to another participant, which has to call requestConnection().
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  void acceptSlavesConnection(const std::string &acceptorName,
                              const std::string &requesterName);

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  void requestSlavesConnection(const std::string &acceptorName,
                               const std::string &requesterName);

    /**
   * Same as acceptSlavesConnection except this only creates the channels, 
   * no vertex list needed!
   */
  void acceptSlavesPreConnection(const std::string &acceptorName,
                                 const std::string &requesterName);

  /** 
   * Same as requestSlavesConnection except this only creates the channels, 
   * no vertex list needed!
   */
  void requestSlavesPreConnection(const std::string &acceptorName,
                                  const std::string &requesterName);

  /*
   * @brief After preliminary communication channels were set up and after 
   *        the mesh partitions were communicated locally for every mesh, 
   *        call this function to update and complete the communication 
   *        channels for every communicated mesh
   */  
  void completeSlavesConnection();

  /**
   * @brief prepares to establish the connections
   *
   * This should be called before calling the accept and request methods.
   * Calling this function forwards the call to the configured master communication.
   *
   * @see com::Communication::prepareEstablishment()
   * @see cleanupEstablishment()
   */
  void prepareEstablishment();

  /**
   * @brief cleans-up to establish the connections
   *
   * This should be called after calling the accept and request methods.
   * Calling this function forwards the call to the configured master communication.
   *
   * @see com::Communication::cleanupEstablishment()
   * @see prepareEstablishment()
   */
  void cleanupEstablishment();

  /**
   * @brief Disconnects from communication space, i.e. participant.
   *
   * This method is called on destruction.
   */
  void closeConnection();

  /// Get the basic communication between the 2 masters.
  com::PtrCommunication getMasterCommunication();

  /// Creates a new distributes communication for that mesh, stores the pointer in _distComs
  void createDistributedCommunication(mesh::PtrMesh mesh);

  /// Sends an array of double values from all slaves (different for each slave).
  void send(double const *itemsToSend,
            int     size,
            int     meshID,
            int     valueDimension);

  /**
   * @brief The master sends a bool to the other master, for performance reasons, we
   * neglect the gathering and checking step.
   */
  void send(bool itemToSend);

  /**
   * @brief The master sends a double to the other master, for performance reasons, we
   * neglect the gathering and checking step.
   */
  void send(double itemToSend);
   
  /// each rank sends its mesh partition to connected ranks  
  void broadcastSendLocalMesh(mesh::Mesh &mesh);

  /// each rank sends the local communication map to the remote connecetd ranks (of the other participant)  
  void broadcastSendLCM(std::map<int, std::vector<int>> &localCommunicationMap, mesh::Mesh &mesh);

  /// each rank sends an int to the remote connected ranks  
  void broadcastSend(int &itemToSend, mesh::Mesh &mesh);

  /// All slaves receive an array of doubles (different for each slave).
  void receive(double *itemsToReceive,
               int     size,
               int     meshID,
               int     valueDimension);

  /// All slaves receive a bool (the same for each slave).
  void receive(bool &itemToReceive);

  /// All slaves receive a double (the same for each slave).
  void receive(double &itemToReceive);

  /// each rank receives mesh partition from connected ranks
  void broadcastReceiveLocalMesh(mesh::Mesh &mesh);

  /// each rank receives local communication maps from remote connetcetd ranks (of the other participant)  
  void broadcastReceiveLCM(std::map<int, std::vector<int>> &localCommunicationMap, mesh::Mesh &mesh);

  /// each rank receives an int from remote connetcetd ranks
  void broadcastReceiveAll(std::vector<int> &itemToReceive, mesh::Mesh &mesh);
  
private:
  logging::Logger _log{"m2n::M2N"};

  /// mesh::getID() -> Pointer to distributed communication
  std::map<int, DistributedCommunication::SharedPointer> _distComs;

  com::PtrCommunication _masterCom;

  DistributedComFactory::SharedPointer _distrFactory;

  bool _isMasterConnected = false;

  bool _areSlavesConnected = false;
};

} // namespace m2n
} // namespace precice
