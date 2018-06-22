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
   * @param[in] nameReuester Name of calling participant.
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
   * @param[in] nameReuester Name of calling participant.
   */
  void requestSlavesConnection(const std::string &acceptorName,
                               const std::string &requesterName);

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
  void send(double *itemsToSend,
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

  /// All slaves receive an array of doubles (different for each slave).
  void receive(double *itemsToReceive,
               int     size,
               int     meshID,
               int     valueDimension);

  /// All slaves receive a bool (the same for each slave).
  void receive(bool &itemToReceive);

  /// All slaves receive a double (the same for each slave).
  void receive(double &itemToReceive);

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
