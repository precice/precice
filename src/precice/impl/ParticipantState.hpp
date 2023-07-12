#pragma once

#include <boost/range/adaptor/map.hpp>
#include <cmath>
#include <memory>
#include <stddef.h>
#include <string>
#include <string_view>
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
#include "precice/impl/ReadGlobalDataContext.hpp"
#include "precice/impl/WriteDataContext.hpp"
#include "precice/impl/WriteGlobalDataContext.hpp"
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

/// Type that represent a compound key of two values
template <typename T>
struct MeshDataKey {
  T mesh;
  T data;
  template <typename Other>
  bool operator<(const MeshDataKey<Other> &other) const
  {
    if (mesh < other.mesh) {
      return true;
    }
    if (other.mesh < mesh) {
      return false;
    }
    return data < other.data;
  }
};

/// Deduction guide for two identical parameter types
template <class T>
MeshDataKey(T, T)->MeshDataKey<T>;

/// Holds coupling state of one participating solver in coupled simulation.
class ParticipantState {
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
  ParticipantState(
      std::string                 name,
      mesh::PtrMeshConfiguration &meshConfig);

  virtual ~ParticipantState();

  /// @name Configuration interface
  /// @{
  /// Adds a configured write \ref Data to the ParticipantState
  void addWriteData(
      const mesh::PtrData &data,
      const mesh::PtrMesh &mesh);

  /// Adds a configured read \ref Data to the ParticipantState
  void addReadData(
      const mesh::PtrData &data,
      const mesh::PtrMesh &mesh);

  /// Adds a configured \ref GlobalData to the ParticipantState
  void addReadGlobalData(
      const mesh::PtrData &data);

  /// Adds a configured write \ref GlobalData to the ParticipantState
  void addWriteGlobalData(
      const mesh::PtrData &data);

  /// Adds a configured read \ref Mapping to the ParticipantState
  void addReadMappingContext(const MappingContext &mappingContext);

  /// Adds a configured write \ref Mapping to the ParticipantState
  void addWriteMappingContext(const MappingContext &mappingContext);

  /// Adds a configured \ref WatchPoint to the ParticipantState
  void addWatchPoint(const PtrWatchPoint &watchPoint);

  /// Adds a configured \ref WatchIntegral to the ParticipantState
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
  void provideMesh(const mesh::PtrMesh &mesh);

  /// Adds a mesh to be received by the participant.
  void receiveMesh(const mesh::PtrMesh &                         mesh,
                   const std::string &                           fromParticipant,
                   double                                        safetyFactor,
                   partition::ReceivedPartition::GeometricFilter geoFilter,
                   const bool                                    allowDirectAccess);
  /// @}

  /// @name Data queries
  /// @{
  /** Provides access to \ref ReadDataContext
   * @pre there exists a \ref ReadDataContext for \ref data
   */
  const ReadDataContext &readDataContext(std::string_view mesh, std::string_view data) const;

  /** Provides access to \ref ReadDataContext
   * @pre there exists a \ref ReadDataContext for \ref data
   */
  ReadDataContext &readDataContext(std::string_view mesh, std::string_view data);

  /**
   * @brief Returns the mesh associated with ReadDataContext with given data name in _readDataContexts of this Participant
   *
   * @param data name of the data
   * @return mesh::PtrMesh, returns nullptr, if no read data contest for given data name was found
   */
  mesh::PtrMesh findMesh(std::string_view data) const;

  /** Provides access to \ref WriteDataContext
   * @pre there exists a \ref WriteDataContext for \ref data
   */
  const WriteDataContext &writeDataContext(std::string_view mesh, std::string_view data) const;

  /** Provides access to \ref WriteDataContext
   * @pre there exists a \ref WriteDataContext for \ref data
   */
  WriteDataContext &writeDataContext(std::string_view mesh, std::string_view data);

  /** Provides access to \ref WriteGlobalDataContext
   * @pre there exists a \ref WriteGlobalDataContext for \ref data
   */
  const WriteGlobalDataContext &writeGlobalDataContext(std::string_view data) const;

  /** Provides access to \ref WriteGlobalDataContext
   * @pre there exists a \ref WriteGlobalDataContext for \ref data
   */
  WriteGlobalDataContext &writeGlobalDataContext(std::string_view data);

  /** Provides access to \ref WriteGlobalDataContext
   * @pre there exists a \ref WriteGlobalDataContext for \ref data
   */
  const ReadGlobalDataContext &readGlobalDataContext(std::string_view data) const;

  /** Provides access to \ref WriteGlobalDataContext
   * @pre there exists a \ref WriteGlobalDataContext for \ref data
   */
  ReadGlobalDataContext &readGlobalDataContext(std::string_view data);

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

  /** Provides const access to all \ref WriteGlobalDataContext objects
   * @remarks does not contain nullptr.
   */
  auto writeGlobalDataContexts() const
  {
    return _writeGlobalDataContexts | boost::adaptors::map_values;
  }

  /** Provides access to all \ref WriteGlobalDataContext objects
   * @remarks does not contain nullptr.
   */
  auto writeGlobalDataContexts()
  {
    return _writeGlobalDataContexts | boost::adaptors::map_values;
  }

  /** Provides const access to all \ref ReadGlobalDataContext objects
   * @remarks does not contain nullptr.
   */
  auto readGlobalDataContexts() const
  {
    return _readGlobalDataContexts | boost::adaptors::map_values;
  }

  /** Provides access to all \ref ReadGlobalDataContext objects
   * @remarks does not contain nullptr.
   */
  auto readGlobalDataContexts()
  {
    return _readGlobalDataContexts | boost::adaptors::map_values;
  }

  /// Is the dataID know to preCICE?
  bool hasData(std::string_view mesh, std::string_view data) const;

  /// Is the data used by this participant?
  bool isDataUsed(std::string_view mesh, std::string_view data) const;

  /// Is the participant allowed to read the data?
  bool isDataRead(std::string_view mesh, std::string_view data) const;

  /// Is the participant allowed to write the data?
  bool isDataWrite(std::string_view mesh, std::string_view data) const;
  /// @}

  /// Is the participant allowed to read the global data?
  bool isGlobalDataRead(std::string_view data) const;

  /// Is the participant allowed to write the global data?
  bool isGlobalDataWrite(std::string_view data) const;

  /// @name Mesh queries
  /// @{
  /*** Provides direct access to a \ref MeshContext given the \ref meshid
   * @param[in] meshID the id of the \ref Mesh
   * @returns a reference to the matching \ref MeshContext
   * @pre the \ref Mesh with \ref meshID is used by the ParticipantState
   */
  const MeshContext &meshContext(std::string_view mesh) const;

  /*** Provides direct access to a \ref MeshContext given the \ref meshid
   * @param[in] meshID the id of the \ref Mesh
   * @returns a reference to the matching \ref MeshContext
   * @pre the \ref Mesh with \ref meshID is used by the ParticipantState
   */
  MeshContext &meshContext(std::string_view mesh);

  /** Provides unordered access to all \ref MeshContext.used by this \ref ParticipantState
   * @remarks The sequence does not contain nullptr
   */
  const std::vector<MeshContext *> &usedMeshContexts() const;

  /** Provides unordered access to all \ref MeshContext.used by this \ref ParticipantState
   * @remarks The sequence does not contain nullptr
   */
  std::vector<MeshContext *> &usedMeshContexts();

  /** Looks for a used MeshContext with a given mesh name.
   * @param[in] name the name of the \ref Mesh
   * @return a reference to the MeshContext
   * @pre there is a matching mesh
   */
  MeshContext &usedMeshContext(std::string_view name);

  /** Looks for a used MeshContext with a given mesh name.
   * @param[in] name the name of the \ref Mesh
   * @return a reference to the MeshContext
   * @pre there is a matching mesh
   */
  MeshContext const &usedMeshContext(std::string_view name) const;

  /// Does preCICE know a mesh with this name?
  bool hasMesh(std::string_view mesh) const;

  /// Is a mesh with this name used by this participant?
  bool isMeshUsed(std::string_view mesh) const;

  /// Is a mesh with this name provided by this participant?
  bool isMeshProvided(std::string_view mesh) const;

  /// Is a mesh with this name received by this participant?
  bool isMeshReceived(std::string_view mesh) const;

  /// Returns whether we are allowed to access a received mesh direct
  /// which requires the config tag <receive-mesh ... direct-access="true"
  bool isDirectAccessAllowed(std::string_view mesh) const;
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
  /// @}

  /// @name Error helpers
  /// @{
  std::string hintForMesh(std::string_view mesh) const;
  std::string hintForMeshData(std::string_view mesh, std::string_view data) const;
  /// @}

private:
  mutable logging::Logger _log{"impl::ParticipantState"};

  std::string _name;

  std::vector<PtrWatchPoint> _watchPoints;

  std::vector<PtrWatchIntegral> _watchIntegrals;

  /// Export contexts to export meshes, data, and more.
  std::vector<io::ExportContext> _exportContexts;

  std::vector<action::PtrAction> _actions;

  template <typename T>
  using MeshMap = std::map<std::string, T, std::less<>>;

  template <typename T>
  using DataMap = std::map<MeshDataKey<std::string>, T, std::less<>>;

  /// All mesh contexts involved in a simulation
  MeshMap<MeshContext *> _meshContexts;

  /// Read mapping contexts used by the participant.
  std::vector<MappingContext> _readMappingContexts;

  /// Write mapping contexts used by the participant.
  std::vector<MappingContext> _writeMappingContexts;

  /// Mesh contexts used by the participant.
  std::vector<MeshContext *> _usedMeshContexts;

  DataMap<WriteDataContext> _writeDataContexts;

  DataMap<ReadDataContext> _readDataContexts;

  std::map<std::string, ReadGlobalDataContext> _readGlobalDataContexts;

  std::map<std::string, WriteGlobalDataContext> _writeGlobalDataContexts;

  bool _useIntraComm = false;

  std::unique_ptr<utils::ManageUniqueIDs> _meshIdManager;

  template <typename ELEMENT_T>
  bool isDataValid(
      const std::vector<ELEMENT_T> &data,
      const ELEMENT_T &             newElement) const;

  void checkDuplicatedUse(std::string_view mesh);

  void checkDuplicatedData(std::string_view mesh, std::string_view data);

  void checkDuplicatedGlobalData(std::string_view data);

  /// To allow white box tests.
  friend struct Integration::Serial::Whitebox::TestConfigurationPeano;
  friend struct Integration::Serial::Whitebox::TestConfigurationComsol;
};

// --------------------------------------------------------- HEADER DEFINITIONS

template <typename ELEMENT_T>
bool ParticipantState::isDataValid(
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
