#pragma once

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "action/Action.hpp"
#include "boost/noncopyable.hpp"
#include "com/Communication.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "m2n/BoundM2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "precice/Participant.hpp"
#include "precice/impl/DataContext.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/Types.hpp"
#include "utils/MultiLock.hpp"

namespace precice {

namespace profiling {
class Event;
}

namespace config {
class Configuration;
}
} // namespace precice

// Forward declaration to friend the boost test struct

namespace Integration::Serial::Whitebox {
struct TestConfigurationPeano;
struct TestConfigurationComsol;
} // namespace Integration::Serial::Whitebox

namespace precice {
namespace cplscheme {
class CouplingSchemeConfiguration;
} // namespace cplscheme
namespace mesh {
class Mesh;
} // namespace mesh

namespace impl {

/// Implementation of Participant. See also pimpl ideom (https://en.cppreference.com/w/cpp/language/pimpl).
class ParticipantImpl {
public:
  ///@name Construction and Configuration
  ///@{

  /**
   * @brief Generic constructor for ParticipantImpl.
   *
   * Use the parameter communicator to specify a custom global MPI communicator.
   * Pass std::nullopt to signal preCICE to use MPI_COMM_WORLD.
   *
   * @param[in] participantName Name of the participant using the interface. Has to
   *        match the name given for a participant in the xml configuration file.
   * @param[in] configurationFileName Name (with path) of the xml configuration file.
   * @param[in] solverProcessIndex If the solver code runs with several processes,
   *        each process using preCICE has to specify its index, which has to start
   *        from 0 and end with solverProcessSize - 1.
   * @param[in] solverProcessSize The number of solver processes using preCICE.
   * @param[in] communicator An optional pointer to an MPI_Comm to use as communicator.
   */
  ParticipantImpl(
      std::string_view      participantName,
      std::string_view      configurationFileName,
      int                   solverProcessIndex,
      int                   solverProcessSize,
      std::optional<void *> communicator);

  /**
   * @brief Destructor
   *
   * Ensures that finalize() has been called.
   *
   * @see finalize
   */
  ~ParticipantImpl();

  ///@}

  /// @name Steering Methods
  ///@{

  /// @copydoc Participant::initialize
  void initialize();

  /// @copydoc Participant::advance
  void advance(double computedTimeStepSize);

  /// @copydoc Participant::finalize
  void finalize();

  ///@}

  ///@name Status Queries
  ///@{

  /// @copydoc Participant::getMeshDimensions
  int getMeshDimensions(std::string_view meshName) const;

  /// @copydoc Participant::getDataDimensions
  int getDataDimensions(std::string_view meshName, std::string_view dataName) const;

  /// @copydoc Participant::isCouplingOngoing
  bool isCouplingOngoing() const;

  /// @copydoc Participant::isTimeWindowComplete
  bool isTimeWindowComplete() const;

  /// @copydoc Participant::getMaxTimeStepSize
  double getMaxTimeStepSize() const;

  ///@}

  ///@name Requirements
  ///@{

  /// @copydoc Participant::requiresInitialData
  bool requiresInitialData();

  /// @copydoc Participant::requiresReadingCheckpoint
  bool requiresReadingCheckpoint();

  /// @copydoc Participant::requiresWritingCheckpoint
  bool requiresWritingCheckpoint();

  ///@}

  ///@name Mesh Access
  ///@anchor precice-mesh-access
  ///@{

  /// @copydoc Participant::resetMesh
  void resetMesh(std::string_view meshName);

  /// @copydoc Participant::resetMeshAccessRegion
  void resetMeshAccessRegion(std::string_view meshName);

  /// @copydoc Participant::requiresMeshConnectivityFor
  bool requiresMeshConnectivityFor(std::string_view meshName) const;

  /// @copydoc Participant::requiresGradientDataFor
  bool requiresGradientDataFor(std::string_view meshName,
                               std::string_view dataName) const;

  /// @copydoc Participant::setMeshVertex
  VertexID setMeshVertex(
      std::string_view              meshName,
      ::precice::span<const double> position);

  /// @copydoc Participant::getMeshVertexSize
  int getMeshVertexSize(std::string_view meshName) const;

  /// @copydoc Participant::setMeshVertices
  void setMeshVertices(
      std::string_view              meshName,
      ::precice::span<const double> positions,
      ::precice::span<VertexID>     ids);

  /// @copydoc Participant::setMeshEdge
  void setMeshEdge(
      std::string_view meshName,
      VertexID         first,
      VertexID         second);

  /// @copydoc Participant::setMeshEdges
  void setMeshEdges(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  /// @copydoc Participant::setMeshTriangle
  void setMeshTriangle(
      std::string_view meshName,
      VertexID         first,
      VertexID         second,
      VertexID         third);

  /// @copydoc Participant::setMeshTriangles
  void setMeshTriangles(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  /// @copydoc Participant::setMeshQuad
  void setMeshQuad(
      std::string_view meshName,
      VertexID         first,
      VertexID         second,
      VertexID         third,
      VertexID         fourth);

  /// @copydoc Participant::setMeshQuads
  void setMeshQuads(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  /// @copydoc Participant::setMeshTetrahedron
  void setMeshTetrahedron(
      std::string_view meshName,
      VertexID         first,
      VertexID         second,
      VertexID         third,
      VertexID         fourth);

  /// @copydoc Participant::setMeshTetrahedra
  void setMeshTetrahedra(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  ///@}

  ///@name Data Access
  ///@{

  /// @copydoc Participant::readData
  void readData(
      std::string_view                meshName,
      std::string_view                dataName,
      ::precice::span<const VertexID> vertices,
      double                          relativeReadTime,
      ::precice::span<double>         values) const;

  /// @copydoc Participant::mapAndReadData
  void mapAndReadData(
      std::string_view              meshName,
      std::string_view              dataName,
      ::precice::span<const double> coordinates,
      double                        relativeReadTime,
      ::precice::span<double>       values) const;

  /// @copydoc Participant::writeAndMapData
  void writeAndMapData(
      std::string_view              meshName,
      std::string_view              dataName,
      ::precice::span<const double> coordinates,
      ::precice::span<const double> values);

  /// @copydoc Participant::writeData
  void writeData(
      std::string_view                meshName,
      std::string_view                dataName,
      ::precice::span<const VertexID> vertices,
      ::precice::span<const double>   values);

  /// @copydoc Participant::writeGradientData
  void writeGradientData(
      std::string_view                meshName,
      std::string_view                dataName,
      ::precice::span<const VertexID> vertices,
      ::precice::span<const double>   gradients);

  ///@}

  /** @name Direct Access
   */
  ///@{

  /// @copydoc Participant::setMeshAccessRegion
  void setMeshAccessRegion(std::string_view              meshName,
                           ::precice::span<const double> boundingBox) const;

  /// @copydoc Participant::getMeshVertexIDsAndCoordinates
  void getMeshVertexIDsAndCoordinates(
      std::string_view          meshName,
      ::precice::span<VertexID> ids,
      ::precice::span<double>   coordinates) const;

  ///@}

  /** @name User profiling
   */
  ///@{

  /// @copydoc Participant::startProfilingSection()
  void startProfilingSection(std::string_view eventName);

  /// @copydoc Participant::stopLastProfilingSection()
  void stopLastProfilingSection();

  ///@}

  /**
   * @brief Allows to access a registered mesh
   */
  /// @todo try to remove or make private. See https://github.com/precice/precice/issues/1269
  const mesh::Mesh &mesh(const std::string &meshName) const;

  struct MappedSamples {
    int write, read;
  };

  /// Returns the amount of mapped read and write samples in the last call to advance.
  MappedSamples mappedSamples() const;

  /// Disable copy construction
  ParticipantImpl(ParticipantImpl const &) = delete;

  /// Disable assignment construction
  ParticipantImpl &operator=(ParticipantImpl const &) = delete;

  /// Disable move construction
  ParticipantImpl(ParticipantImpl &&) = delete;

  /// Disable move assignment
  ParticipantImpl &operator=(ParticipantImpl &&) = delete;

private:
  mutable logging::Logger _log{"impl::ParticipantImpl"};

  std::string _accessorName;

  int _accessorProcessRank;

  int _accessorCommunicatorSize;

  impl::PtrParticipant _accessor;

  /// Spatial dimensions of problem.
  int _dimensions = 0;

  utils::MultiLock<std::string> _meshLock;

  /// mesh name to mesh ID mapping.
  std::map<std::string, int> _meshIDs;

  std::map<std::string, m2n::BoundM2N> _m2ns;

  /// Holds information about solvers participating in the coupled simulation.
  std::vector<impl::PtrParticipant> _participants;

  cplscheme::PtrCouplingScheme _couplingScheme;

  /// Represents the various states a Participant can be in.
  enum struct State {
    Constructed, // Initial state of Participant
    Initialized, // Participant.initialize() triggers transition from State::Constructed to State::Initialized; mandatory
    Finalized    // Participant.finalize() triggers transition form State::Initialized to State::Finalized; mandatory
  };

  /// Are experimental API calls allowed?
  bool _allowsExperimental = false;

  /// Are experimental remeshing API calls allowed?
  bool _allowsRemeshing = false;

  /// Are participants waiting for each other in finalize?
  bool _waitInFinalize = false;

  /// The current State of the Participant
  State _state{State::Constructed};

  /// Counts calls to advance for plotting.
  long int _numberAdvanceCalls = 0;

  /// Counts the amount of samples mapped in write mappings executed in the latest advance
  int _executedWriteMappings = 0;

  /// Counts the amount of samples mapped in read mappings executed in the latest advance
  int _executedReadMappings = 0;

  /// The hash of the configuration file used to configure this participant
  std::string _configHash;

  /**
   * @brief Configures the coupling interface from the given xml file.
   *
   * Only after the configuration a reasonable state of a ParticipantImpl
   * object is achieved.
   *
   * @param[in] configurationFileName Name (with path) of the xml config. file.
   */
  void configure(std::string_view configurationFileName);

  /**
   * @brief Configures the coupling interface with a prepared configuration.
   *
   * Can be used to configure the ParticipantImpl without xml file. Requires
   * to manually setup the configuration object.
   */
  void configure(const config::Configuration &configuration);

  void configureM2Ns(const m2n::M2NConfiguration::SharedPointer &config);

  enum struct ExportTiming : bool {
    Advance = false,
    Initial = true
  };

  /// Exports meshes with data and watch point data.
  /// @param[in] timing when the exports are requested
  void handleExports(ExportTiming timing);

  /// Determines participants providing meshes to other participants.
  void configurePartitions(
      const m2n::M2NConfiguration::SharedPointer &m2nConfig);

  /// Communicate bounding boxes and look for overlaps
  void compareBoundingBoxes();

  /// Communicate meshes and create partitions
  void computePartitions();

  /// Helper for mapWrittenData and mapReadData
  void computeMappings(std::vector<MappingContext> &contexts, const std::string &mappingType);

  /// Computes, and performs write mappings of the initial data in initialize
  void mapInitialWrittenData();

  /// Computes, and performs suitable write mappings either entirely or after given time
  void mapWrittenData(std::optional<double> after = std::nullopt);

  // Computes, and performs read mappings of the initial data in initialize
  void mapInitialReadData();

  // Computes, and performs read mappings
  void mapReadData();

  /**
   * @brief Removes samples in mapped to data connected to received data via a mapping.
   *
   * This prevents old samples from blocking remappings.
   */
  void trimReadMappedData(double timeAfterAdvance, bool isTimeWindowComplete, const cplscheme::ImplicitData &fromData);

  /**
   * @brief Performs all data actions with given timing.
   *
   * @param[in] timings the timings of the action.
   */
  void performDataActions(const std::set<action::Action::Timing> &timings);

  /**
   * @brief Resets written data.
   *
   * @param isAtWindowEnd set true, if function is called at end of window to also trim the time sample storage
   * @param isTimeWindowComplete set true, if function is called at end of converged window to trim and move the sample storage.
   */
  void resetWrittenData(); //bool isAtWindowEnd, bool isTimeWindowComplete);

  /// Determines participant accessing this interface from the configuration.
  impl::PtrParticipant determineAccessingParticipant(
      const config::Configuration &config);

  /// Initializes intra-participant communication.
  void initializeIntraCommunication();

  /// Advances the coupling schemes
  void advanceCouplingScheme();

  /// Syncs the time step size between all ranks (all time steps sizes should be the same!)
  void syncTimestep(double computedTimeStepSize);

  /// Which channels to close in closeCommunicationChannels()
  enum class CloseChannels : bool {
    All         = false,
    Distributed = true
  };

  /// Syncs the primary ranks of all connected participants
  void closeCommunicationChannels(CloseChannels cc);

  /// Completes everything data-related between adding time to and advancing the coupling scheme
  void handleDataBeforeAdvance(bool reachedTimeWindowEnd, double timeSteppedTo, bool performedReinit);

  /// Completes everything data-related after advancing the coupling scheme
  void handleDataAfterAdvance(bool reachedTimeWindowEnd, bool isTimeWindowComplete, double timeSteppedTo, double timeAfterAdvance, const cplscheme::ImplicitData &receivedData);

  /// Creates a Stample at the given time for each write Data and zeros the buffers
  void samplizeWriteData(double time, bool performedReinit = false);

  /// Discards data before the given time for all meshes and data known by this participant
  void trimOldDataBefore(double time);

  /// Discards send (currently write) data of a participant after a given time when another iteration is required
  void trimSendDataAfter(double time);

  /// How many ranks have changed each used mesh
  using MeshChanges = std::vector<int>;

  /** Allreduce of the amount of changed meshes on each rank.
   * @return a vector of the size of meshcontexts which contain the amount of ranks that changed each mesh
   */
  MeshChanges getTotalMeshChanges() const;

  /// Clears stample of changed meshes to make them consistent after the reinitialization.
  void clearStamplesOfChangedMeshes(MeshChanges totalMeshChanges);

  /** Exchanges request to remesh with all connecting participants.
   *
   * @param[in] requestReinit does this participant request to remesh?
   *
   * @return does any participant request to remesh?
   */
  bool reinitHandshake(bool requestReinit) const;

  /// Reinitializes preCICE
  void reinitialize();

  /// Connect participants including repartitioning
  void setupCommunication();

  /// Setup mesh watcher such as WatchPoints
  void setupWatcher();

  /// Returns if a user has to define an access region for direct
  /// mesh access and just-in-time mapping or not
  /// Right now, that's required in parallel runs on received meshes
  bool requiresUserDefinedAccessRegion(std::string_view meshName) const;

  /// To allow white box tests.
  friend struct Integration::Serial::Whitebox::TestConfigurationPeano;
  friend struct Integration::Serial::Whitebox::TestConfigurationComsol;

  std::unique_ptr<profiling::Event> _solverInitEvent;
  std::unique_ptr<profiling::Event> _solverAdvanceEvent;

  std::vector<profiling::Event> _userEvents;
};

} // namespace impl
} // namespace precice
