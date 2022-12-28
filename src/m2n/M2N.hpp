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

// Forward declaration to friend unit tests which only use the primary com
struct WhiteboxAccessor;

/**
 * @brief M2N communication class.
 * This layer is necessary since communication between two participants can be working via several meshes,
 * each possibly with a different decomposition. In principle, this class is only a map from meshes to DistributedCommunications
 *
 */
class M2N {
public:
  M2N(com::PtrCommunication intraComm, DistributedComFactory::SharedPointer distrFactory, bool useOnlyPrimaryCom = false, bool useTwoLevelInit = false);

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
  void acceptPrimaryRankConnection(const std::string &acceptorName,
                                   const std::string &requesterName);

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  void requestPrimaryRankConnection(const std::string &acceptorName,
                                    const std::string &requesterName);

  /**
   * @brief Connects to another participant, which has to call requestConnection().
   *
   * @param[in] acceptorName Name of calling participant.
   * @param[in] requesterName Name of remote participant to connect to.
   */
  void acceptSecondaryRanksConnection(const std::string &acceptorName,
                                      const std::string &requesterName);

  /**
   * @brief Connects to another participant, which has to call acceptConnection().
   *
   * @param[in] acceptorName Name of remote participant to connect to.
   * @param[in] requesterName Name of calling participant.
   */
  void requestSecondaryRanksConnection(const std::string &acceptorName,
                                       const std::string &requesterName);

  /**
   * Same as acceptSecondaryRanksConnection except this only creates the channels,
   * no vertex list needed!
   */
  void acceptSecondaryRanksPreConnection(const std::string &acceptorName,
                                         const std::string &requesterName);

  /**
   * Same as requestSecondaryRanksConnection except this only creates the channels,
   * no vertex list needed!
   */
  void requestSecondaryRanksPreConnection(const std::string &acceptorName,
                                          const std::string &requesterName);

  /*
   * @brief After preliminary communication channels were set up and after
   *        the mesh partitions were communicated locally for every mesh,
   *        call this function to update and complete the communication
   *        channels for every communicated mesh
   */
  void completeSecondaryRanksConnection();

  /**
   * @brief prepares to establish the connections
   *
   * This should be called before calling the accept and request methods.
   * Calling this function forwards the call to the configured primary communication.
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
   * Calling this function forwards the call to the configured primary communication.
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
   * Calls closePrimaryRankConnection() and closeSecondaryRanksConnection()
   * This method is called on destruction.
   */
  void closeConnection();

  /// Disconnects the primary connection
  void closePrimaryRankConnection();

  /// Disconnects all connections of the DistributedCommunication
  void closeDistributedConnections();

  /// Get the basic communication between the 2 primary ranks.
  com::PtrCommunication getPrimaryRankCommunication();

  /// Creates a new distributes communication for that mesh, stores the pointer in _distComs
  void createDistributedCommunication(const mesh::PtrMesh &mesh);

  /// Sends an array of double values from all ranks (different for each rank).
  void send(precice::span<double const> itemsToSend,
            int                         meshID,
            int                         valueDimension);

  /**
   * @brief The primary rank sends a bool to the other primary rank, for performance reasons, we
   * neglect the gathering and checking step.
   */
  void send(bool itemToSend);

  /**
   * @brief The primary rank sends a double to the other primary rank, for performance reasons, we
   * neglect the gathering and checking step.
   */
  void send(double itemToSend);

  /**
   * @brief The primary rank sends an int to the other primary rank, for performance reasons, we
   * neglect the gathering and checking step.
   */
  void send(int itemToSend);

  /// Broadcasts a mesh to connected ranks on remote participant (concerning the given mesh)
  void broadcastSendMesh(mesh::Mesh &mesh);

  /// Scatters a communication map over connected ranks on remote participant (concerning the given mesh)
  void scatterAllCommunicationMap(std::map<int, std::vector<int>> &localCommunicationMap, mesh::Mesh &mesh);

  /// Broadcasts an int to connected ranks on remote participant (concerning the given mesh)
  void broadcastSend(int itemToSend, mesh::Mesh &mesh);

  /// All ranks receive an array of doubles (different for each rank).
  /// The values received can be gradient data
  /// Gradient dimension : 0: dx-values, 1: dy-values, 2:dz-values
  void receive(precice::span<double> itemsToReceive,
               int                   meshID,
               int                   valueDimension);

  /// All ranks receive a bool (the same for each rank).
  void receive(bool &itemToReceive);

  /// All ranks receive a double (the same for each rank).
  void receive(double &itemToReceive);

  /// All ranks receive an int (the same for each rank).
  void receive(int &itemToReceive);

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

  com::PtrCommunication _intraComm;

  DistributedComFactory::SharedPointer _distrFactory;

  bool _isPrimaryRankConnected = false;

  bool _areSecondaryRanksConnected = false;

  // The following flag is (solely) needed for unit tests between two serial participants.
  // To also use the secondary communication would require a lengthy setup of meshes
  // and their re-partitioning, which could also not be moved to some fixture as the M2Ns
  // are created through the configuration.
  // See e.g. "CplSchemeTests/ExplicitCouplingSchemeTests/testConfiguredSimpleExplicitCoupling"
  // This flag gives a loophole. It is set to false for normal use and modified in the
  // respective tests through a friend declaration.

  /// between two serial participants, only use the primary com and no secondary com
  bool _useOnlyPrimaryCom = false;

  /// use the two-level initialization concept
  bool _useTwoLevelInit = false;

  // @brief To allow access to _useOnlyPrimaryCom
  friend struct WhiteboxAccessor;
};

/// struct giving access _useOnlyPrimaryCom
struct WhiteboxAccessor {
  static auto useOnlyPrimaryCom(PtrM2N m2n) -> typename std::add_lvalue_reference<decltype(m2n->_useOnlyPrimaryCom)>::type
  {
    return m2n->_useOnlyPrimaryCom;
  }
};

} // namespace m2n
} // namespace precice
