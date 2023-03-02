#pragma once

#include <Eigen/Core>
#include <boost/range/adaptor/map.hpp>
#include <cmath>
#include <memory>
#include <stddef.h>
#include <string>
#include <utility>
#include <vector>

#include "SharedPointer.hpp"
#include "action/SharedPointer.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "io/ExportContext.hpp"
#include "io/config/ExportConfiguration.hpp"
#include "logging/Logger.hpp"
#include "mapping/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "partition/ReceivedPartition.hpp"
#include "precice/impl/ReadDataContext.hpp"
#include "precice/impl/WriteDataContext.hpp"
#include "precice/types.hpp"
#include "utils/IntraComm.hpp"
#include "utils/ManageUniqueIDs.hpp"

namespace precice {
namespace impl {
struct MeshContext;
struct MappingContext;
} // namespace impl
} // namespace precice

// Forward declaration to friend the boost test struct
namespace Integration {
namespace Serial {
namespace Whitebox {
struct TestConfigurationPeano;
struct TestConfigurationComsol;
} // namespace Whitebox
} // namespace Serial
} // namespace Integration

namespace precice {
namespace utils {
class ManageUniqueIDs;
} // namespace utils

namespace impl {

/// Holds coupling state of one participating solver in coupled simulation.
class Participant {
public:
  enum MappingConstants {
    MAPPING_LINEAR_CONSERVATIVE,
    MAPPING_LINEAR_CONSISTENT,
    MAPPING_DIRECT
  };

  /**
   * @brief Constructor.
   *
   * @param[in] name Name of the participant. Has to be unique.
   */
  Participant(
      std::string                 name,
      mesh::PtrMeshConfiguration &meshConfig);

  virtual ~Participant();

  /// @name Configuration interface
  /// @{
  /// Adds a configured write \ref Data to the Participant
  void addWriteData(
      const mesh::PtrData &data,
      const mesh::PtrMesh &mesh);

  /// Adds a configured read \ref Data to the Participant
  void addReadData(
      const mesh::PtrData &data,
      const mesh::PtrMesh &mesh,
      int                  interpolationOrder);

  /// Adds a configured read \ref Mapping to the Participant
  void addReadMappingContext(const MappingContext &mappingContext);

  /// Adds a configured write \ref Mapping to the Participant
  void addWriteMappingContext(const MappingContext &mappingContext);

  /// Adds a configured \ref WatchPoint to the Participant
  void addWatchPoint(const PtrWatchPoint &watchPoint);

  /// Adds a configured \ref WatchIntegral to the Participant
  void addWatchIntegral(const PtrWatchIntegral &watchIntegral);

  /// Sets weather the participant was configured with a primary tag
  void setUsePrimaryRank(bool useIntraComm);

  /// Sets the manager responsible for providing unique IDs to meshes.
  void setMeshIdManager(std::unique_ptr<utils::ManageUniqueIDs> &&idm)
  {
    _meshIdManager = std::move(idm);
  }

  /// Adds a configured \ref Action to the participant
  void addAction(action::PtrAction &&action);

  /// Adds a configured \ref ExportContext to export meshes and data.
  void addExportContext(const io::ExportContext &context);

  /// Adds a mesh to be provided by the participant.
  void provideMesh(const mesh::PtrMesh &mesh, bool dynamic = false);

  /// Adds a mesh to be received by the participant.
  void receiveMesh(const mesh::PtrMesh &                         mesh,
                   const std::string &                           fromParticipant,
                   double                                        safetyFactor,
                   partition::ReceivedPartition::GeometricFilter geoFilter,
                   const bool                                    allowDirectAccess);

  void registerDynamicParticipant(const std::string &name);
  /// @}

  /// @name Data queries
  /// @{
  /** Provides access to \ref ReadDataContext
   * @pre there exists a \ref ReadDataContext for \ref dataID
   */
  const ReadDataContext &readDataContext(DataID dataID) const;

  /** Provides access to \ref ReadDataContext
   * @pre there exists a \ref ReadDataContext for \ref dataID
   */
  ReadDataContext &readDataContext(DataID dataID);

  /**
   * Provides access to \ref ReadDataContext
   * @pre there exists a \ref ReadDataContext for \ref dataName
   */
  ReadDataContext &readDataContext(const std::string &dataName);

  /** Provides access to \ref WriteDataContext
   * @pre there exists a \ref WriteDataContext for \ref dataID
   */
  const WriteDataContext &writeDataContext(DataID dataID) const;

  /** Provides access to \ref WriteDataContext
   * @pre there exists a \ref WriteDataContext for \ref dataID
   */
  WriteDataContext &writeDataContext(DataID dataID);

  /** Provides access to all \ref WriteDataContext objects
   * @remarks does not contain nullptr.
   */
  auto writeDataContexts()
  {
    return _writeDataContexts | boost::adaptors::map_values;
  }

  /** Provides access to all \ref ReadDataContext objects
   * @remarks does not contain nullptr.
   */
  auto readDataContexts()
  {
    return _readDataContexts | boost::adaptors::map_values;
  }

  /** @brief Determines and returns the maximum order of all read waveforms of this participant
   */
  int maxReadWaveformOrder() const
  {
    int maxOrder = -1;
    for (auto &context : _readDataContexts | boost::adaptors::map_values) {
      maxOrder = std::max(maxOrder, context.getInterpolationOrder());
    }
    return maxOrder;
  }

  /// Is the dataID know to preCICE?
  bool hasData(DataID dataID) const;

  /// Is the data used by this participant?
  bool isDataUsed(const std::string &dataName, MeshID meshId) const;

  /// Is the participant allowed to read the data?
  bool isDataRead(DataID dataID) const;

  /// Is the participant allowed to write the data?
  bool isDataWrite(DataID dataID) const;

  /// What is the dataID of the used data from a used mesh given the meshid and the data name?
  int getUsedDataID(const std::string &dataName, MeshID meshID) const;

  /// What is the name of the given data id
  std::string getDataName(DataID dataID) const;
  /// @}

  /// @name Mesh queries
  /// @{
  /*** Provides direct access to a \ref MeshContext given the \ref meshid
   * @param[in] meshID the id of the \ref Mesh
   * @returns a reference to the matching \ref MeshContext
   * @pre the \ref Mesh with \ref meshID is used by the Participant
   */
  const MeshContext &meshContext(MeshID meshID) const;

  /*** Provides direct access to a \ref MeshContext given the \ref meshid
   * @param[in] meshID the id of the \ref Mesh
   * @returns a reference to the matching \ref MeshContext
   * @pre the \ref Mesh with \ref meshID is used by the Participant
   */
  MeshContext &meshContext(MeshID meshID);

  /** Provides unordered access to all \ref MeshContext.used by this \ref Participant
   * @remarks The sequence does not contain nullptr
   */
  const std::vector<MeshContext *> &usedMeshContexts() const;

  /** Provides unordered access to all \ref MeshContext.used by this \ref Participant
   * @remarks The sequence does not contain nullptr
   */
  std::vector<MeshContext *> &usedMeshContexts();

  /** Looks for a used MeshContext with a given mesh name.
   * @param[in] name the name of the \ref Mesh
   * @return a reference to the MeshContext
   * @pre there is a matching mesh
   */
  MeshContext &usedMeshContext(const std::string &name);

  /** Looks for a used MeshContext with a given mesh name.
   * @param[in] name the name of the \ref Mesh
   * @return a reference to the MeshContext
   * @pre there is a matching mesh
   */
  MeshContext const &usedMeshContext(const std::string &name) const;

  /** Looks for a used MeshContext with a given mesh ID.
   * @param[in] meshID the id of the \ref Mesh
   * @return a reference to the MeshContext
   * @pre there is a matching mesh
   */
  MeshContext &usedMeshContext(MeshID meshID);

  /** Looks for a used MeshContext with a given meshID
   * @param[in] meshID the id of the \ref Mesh
   * @return a reference to the MeshContext
   * @pre there is a matching mesh
   */
  MeshContext const &usedMeshContext(MeshID meshID) const;

  /// Does preCICE know a mesh with this meshID?
  bool hasMesh(MeshID meshID) const;

  /// Does preCICE know a mesh with this name?
  bool hasMesh(const std::string &meshName) const;

  /// Is a mesh with this id used by this participant?
  bool isMeshUsed(MeshID meshID) const;

  /// Is a mesh with this name used by this participant?
  bool isMeshUsed(const std::string &meshID) const;

  /// Is a mesh with this id provided?
  bool isMeshProvided(MeshID meshID) const;

  /// Is a mesh with this name provided by this participant?
  bool isMeshProvided(const std::string &meshName) const;

  /// Is a mesh with this name received by this participant?
  bool isMeshReceived(const std::string &meshName) const;

  /// Get the used mesh id of a mesh with this name.
  int getUsedMeshID(const std::string &meshName) const;

  /// Returns whether we are allowed to access a received mesh direct
  /// which requires the config tag <receive-mesh ... direct-access="true"
  bool isDirectAccessAllowed(const int meshID) const;

  /// Get the name of a mesh given by its id.
  std::string getMeshName(MeshID meshID) const;

  /// Get a mesh name which uses the given data id.
  std::string getMeshNameFromData(DataID dataID) const;
  /// @}

  /// @name Exporting interface
  /// @{
  /// Exports the initial state of meshes
  void exportInitial();

  /// Exports the final state of meshes
  void exportFinal();

  struct IntermediateExport {
    size_t timewindow;
    size_t iteration;
    double time;
    bool   complete;
  };

  /// Exports timewindows and iterations of meshes and watchpoints
  void exportIntermediate(IntermediateExport exp);

  /// @}

  /// @name Other queries
  /// @{
  /// Returns the name of the participant.
  const std::string &getName() const;

  /// Returns true, if the participant uses a primary tag.
  bool useIntraComm() const;

  /// Provided access to all read \ref MappingContext
  std::vector<MappingContext> &readMappingContexts();

  /// Provided access to all write \ref MappingContext
  std::vector<MappingContext> &writeMappingContexts();

  /// Provided access to all \ref WatchPoints
  std::vector<PtrWatchPoint> &watchPoints();

  /// Provided access to all \ref WatchIntegrals
  std::vector<PtrWatchIntegral> &watchIntegrals();

  /// Provided access to all \ref Action
  std::vector<action::PtrAction> &actions();

  /// Provided access to all \ref Action
  const std::vector<action::PtrAction> &actions() const;

  /// Returns all \ref ExportContext for exporting meshes and data.
  const std::vector<io::ExportContext> &exportContexts() const;

  /// Is this participant dynamic?
  bool isDynamic() const;

  /// Returns the names of all dynamic participants including the local.
  std::set<std::string> dynamicParticipants() const;
  /// @}

private:
  mutable logging::Logger _log{"impl::Participant"};

  std::string _name;

  std::vector<PtrWatchPoint> _watchPoints;

  std::vector<PtrWatchIntegral> _watchIntegrals;

  /// Export contexts to export meshes, data, and more.
  std::vector<io::ExportContext> _exportContexts;

  std::vector<action::PtrAction> _actions;

  /// All mesh contexts involved in a simulation, mesh ID == index.
  std::vector<MeshContext *> _meshContexts; // @todo use map here!

  /// Read mapping contexts used by the participant.
  std::vector<MappingContext> _readMappingContexts;

  /// Write mapping contexts used by the participant.
  std::vector<MappingContext> _writeMappingContexts;

  /// Mesh contexts used by the participant.
  std::vector<MeshContext *> _usedMeshContexts;

  std::map<DataID, WriteDataContext> _writeDataContexts;

  std::map<DataID, ReadDataContext> _readDataContexts;

  bool _useIntraComm = false;

  std::unique_ptr<utils::ManageUniqueIDs> _meshIdManager;

  std::set<std::string> _dynamicParticipants;

  template <typename ELEMENT_T>
  bool isDataValid(
      const std::vector<ELEMENT_T> &data,
      const ELEMENT_T &             newElement) const;

  void checkDuplicatedUse(const mesh::PtrMesh &mesh);

  void checkDuplicatedData(const mesh::PtrData &data, const std::string &meshName);

  /// To allow white box tests.
  friend struct Integration::Serial::Whitebox::TestConfigurationPeano;
  friend struct Integration::Serial::Whitebox::TestConfigurationComsol;
};

// --------------------------------------------------------- HEADER DEFINITIONS

template <typename ELEMENT_T>
bool Participant::isDataValid(
    const std::vector<ELEMENT_T> &data,
    const ELEMENT_T &             newElement) const
{
  for (size_t i = 0; i < data.size(); i++) {
    if (data[i].name == newElement.name) {
      return false;
    }
  }
  return true;
}

} // namespace impl
} // namespace precice
