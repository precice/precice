#pragma once

#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include "DistributedComFactory.hpp"
#include "SharedPointer.hpp"
#include "com/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "m2n/DistributedCommunication.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace mesh {
class Mesh;
} // namespace mesh

namespace m2n {

// Forward declaration to friend unit tests which only use the master com
struct WhiteboxAccessor;

/**
 * @brief M2N communication class.
 * This layer is necessary since communication between two participants can be working via several meshes,
 * each possibly with a different decomposition. In principle, this class is only a map from meshes to DistributedCommunications
 *
 */
class M2N {
public:
  M2N(com::PtrCommunication masterCom, DistributedComFactory::SharedPointer distrFactory, bool useOnlyMasterCom = false, bool useTwoLevelInit = false);

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
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   *
   * @see com::Communication::prepareEstablishment()
   * @see cleanupEstablishment()
   */
  void prepareEstablishment(const std::string &acceptorName,
                            const std::string &requesterName);

  /**
   * @brief cleans-up to establish the connections
   *
   * This should be called after calling the accept and request methods.
   * Calling this function forwards the call to the configured master communication.
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   *
   * @see com::Communication::cleanupEstablishment()
   * @see prepareEstablishment()
   */
  void cleanupEstablishment(const std::string &acceptorName,
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
  void send(double const *itemsToSend,
            int           size,
            int           meshID,
            int           valueDimension);

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

  /// Broadcasts a mesh to connected ranks on remote participant (concerning the given mesh)
  void broadcastSendMesh(mesh::Mesh &mesh);

  /// Scatters a communication map over connected ranks on remote participant (concerning the given mesh)
  void scatterAllCommunicationMap(std::map<int, std::vector<int>> &localCommunicationMap, mesh::Mesh &mesh);

  /// Broadcasts an int to connected ranks on remote participant (concerning the given mesh)
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

  /// Receive mesh partitions per connected rank on remote participant (concerning the given mesh)
  void broadcastReceiveAllMesh(mesh::Mesh &mesh);

  /// Gathers a communication maps from connected ranks on remote participant (concerning the given mesh)
  void gatherAllCommunicationMap(std::map<int, std::vector<int>> &localCommunicationMap, mesh::Mesh &mesh);

  /**
   * @brief Receives an int per connected rank on remote participant (concerning the given mesh)
   * @para[out] itemToReceive received ints from remote ranks are stored with the sender rank order
   */
  void broadcastReceiveAll(std::vector<int> &itemToReceive, mesh::Mesh &mesh);

  bool usesTwoLevelInitialization()
  {
    return _useTwoLevelInit;
  }

private:
  logging::Logger _log{"m2n::M2N"};

  /// mesh::getID() -> Pointer to distributed communication
  std::map<int, DistributedCommunication::SharedPointer> _distComs;

  com::PtrCommunication _masterCom;

  DistributedComFactory::SharedPointer _distrFactory;

  bool _isMasterConnected = false;

  bool _areSlavesConnected = false;

  // The following flag is (solely) needed for unit tests between two serial participants.
  // To also use the slaves-slaves communication would require a lengthy setup of meshes
  // and their re-partitioning, which could also not be moved to some fixture as the M2Ns
  // are created through the configuration.
  // See e.g. "CplSchemeTests/ExplicitCouplingSchemeTests/testConfiguredSimpleExplicitCoupling"
  // This flag gives a loophole. It is set to false for normal use and modfied in the
  // respective tests through a fried decleration.

  /// between two serial participants, only use the master-master com and no slaves-slaves com
  bool _useOnlyMasterCom = false;

  /// use the two-level initialization concept
  bool _useTwoLevelInit = false;

  // @brief To allow access to _useOnlyMasterCom
  friend struct WhiteboxAccessor;
};

/// struct giving access _useOnlyMasterCom
struct WhiteboxAccessor {
  static auto useOnlyMasterCom(PtrM2N m2n) -> typename std::add_lvalue_reference<decltype(m2n->_useOnlyMasterCom)>::type
  {
    return m2n->_useOnlyMasterCom;
  }
};

} // namespace m2n
} // namespace precice
