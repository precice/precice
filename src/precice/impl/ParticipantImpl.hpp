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
#include "cplscheme/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "m2n/BoundM2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "precice/Participant.hpp"
#include "precice/impl/DataContext.hpp"
#include "precice/impl/GlobalDataContext.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/types.hpp"
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
namespace Integration {
namespace Serial {
namespace Whitebox {
struct TestConfigurationPeano;
struct TestConfigurationComsol;
} // namespace Whitebox
} // namespace Serial
} // namespace Integration

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

  /// @copydoc Participant::hasMesh
  bool hasMesh(std::string_view meshName) const;

  /// @copydoc Participant::requiresMeshConnectivityFor
  bool requiresMeshConnectivityFor(std::string_view meshName) const;

  /// @copydoc Participant::requiresGradientDataFor
  bool requiresGradientDataFor(std::string_view meshName,
                               std::string_view dataName) const;

  /// @copydoc Participant::setMeshVertex
  int setMeshVertex(
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
      int              firstVertexID,
      int              secondVertexID);

  /// @copydoc Participant::setMeshEdges
  void setMeshEdges(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  /// @copydoc Participant::setMeshTriangle
  void setMeshTriangle(
      std::string_view meshName,
      int              firstVertexID,
      int              secondVertexID,
      int              thirdVertexID);

  /// @copydoc Participant::setMeshTriangles
  void setMeshTriangles(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  /// @copydoc Participant::setMeshQuad
  void setMeshQuad(
      std::string_view meshName,
      int              firstVertexID,
      int              secondVertexID,
      int              thirdVertexID,
      int              fourthVertexID);

  /// @copydoc Participant::setMeshQuads
  void setMeshQuads(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  /// @copydoc Participant::setMeshTetrahedron
  void setMeshTetrahedron(
      std::string_view meshName,
      int              firstVertexID,
      int              secondVertexID,
      int              thirdVertexID,
      int              fourthVertexID);

  /// @copydoc Participant::setMeshTetrahedra
  void setMeshTetrahedra(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  ///@}

  ///@name Data Access
  ///@{

  /// @copydoc Participant::hasData
  bool hasData(
      std::string_view meshName,
      std::string_view dataName) const;

  /// @copydoc SolverInterface::hasGlobalData
  bool hasGlobalData(const std::string &dataName) const;
  
  /// @copydoc SolverInterface::getGlobalDataID
  int getGlobalDataID(const std::string &dataName) const;
  
  /// @copydoc Participant::readData
  void readData(
      std::string_view                meshName,
      std::string_view                dataName,
      ::precice::span<const VertexID> vertices,
      double                          relativeReadTime,
      ::precice::span<double>         values) const;

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

/// @copydoc SolverInterface::writeGlobalVectorData
  void writeGlobalVectorData(
      int           fromDataID,
      const double *value);


    /// @copydoc SolverInterface::writeGlobalScalarData
  void writeGlobalScalarData(
      int    dataID,
      double value);

  /// @copydoc SolverInterface::readGlobalScalarData(int, double&) const
  void readGlobalScalarData(
      int     toDataID,
      double &value) const;

  /// @copydoc SolverInterface::readGlobalScalarData(int, double, double&) const
  void readGlobalScalarData(
      int     toDataID,
      double  relativeReadTime,
      double &value) const;
  ///@}

  /** @name Experimental Data Access
   * These API functions are \b experimental and may change in future versions.
   */
  ///@{

  /// @copydoc Participant::setMeshAccessRegion
  void setMeshAccessRegion(std::string_view              meshName,
                           ::precice::span<const double> boundingBox) const;

  /// @copydoc Participant::getMeshVerticesAndIDs
  void getMeshVerticesAndIDs(
      std::string_view          meshName,
      ::precice::span<VertexID> ids,
      ::precice::span<double>   coordinates) const;

  ///@}

  /**
   * @brief Allows to access a registered mesh
   */
  /// @todo try to remove or make private. See https://github.com/precice/precice/issues/1269
  const mesh::Mesh &mesh(const std::string &meshName) const;

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

  /// setMeshAccessRegion may only be called once
  mutable bool _accessRegionDefined = false;

  /// The current State of the Participant
  State _state{State::Constructed};

  /// Counts calls to advance for plotting.
  long int _numberAdvanceCalls = 0;

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

  /// Exports meshes with data and watch point data.
  void handleExports();

  /// Determines participants providing meshes to other participants.
  void configurePartitions(
      const m2n::M2NConfiguration::SharedPointer &m2nConfig);

  /// Communicate bounding boxes and look for overlaps
  void compareBoundingBoxes();

  /// Communicate meshes and create partitions
  void computePartitions();

  /// Helper for mapWrittenData and mapReadData
  void computeMappings(std::vector<MappingContext> &contexts, const std::string &mappingType);

  /// Computes, performs, and resets all suitable write mappings.
  void mapWrittenData();

  /// Computes, performs, and resets all suitable read mappings.
  void mapReadData();

  /**
   * @brief Performs all data actions with given timing.
   *
   * @param[in] timings the timings of the action.
   * @param[in] time the current total simulation time.
   */
  void performDataActions(
      const std::set<action::Action::Timing> &timings,
      double                                  time);

  /**
   * @brief Resets written data.
   *
   * @param isAtWindowEnd set true, if function is called at end of window to also trim the time sample storage
   * @param isTimeWindowComplete set true, if function is called at end of converged window to trim and move the sample storage.
   */
  void resetWrittenData(bool isAtWindowEnd, bool isTimeWindowComplete);

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

  void readGlobalVectorDataImpl(
      int     toDataID,
      double  relativeReadTime,
      double *value) const;


  void readGlobalScalarDataImpl(
      int     toDataID,
      double  relativeReadTime,
      double &value) const;

  /// To allow white box tests.
  friend struct Integration::Serial::Whitebox::TestConfigurationPeano;
  friend struct Integration::Serial::Whitebox::TestConfigurationComsol;

  std::unique_ptr<profiling::Event> _solverInitEvent;
  std::unique_ptr<profiling::Event> _solverAdvanceEvent;
};

} // namespace impl
} // namespace precice
