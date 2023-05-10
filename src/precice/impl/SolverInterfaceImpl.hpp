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
#include "precice/SolverInterface.hpp"
#include "precice/impl/DataContext.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/types.hpp"
#include "utils/MultiLock.hpp"

namespace precice {
namespace config {
class SolverInterfaceConfiguration;
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

/// Implementation of SolverInterface. See also pimpl ideom (https://en.cppreference.com/w/cpp/language/pimpl).
class SolverInterfaceImpl {
public:
  ///@name Construction and Configuration
  ///@{

  /**
   * @brief Generic constructor for SolverInterfaceImpl.
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
  SolverInterfaceImpl(
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
  ~SolverInterfaceImpl();

  ///@}

  /// @name Steering Methods
  ///@{

  /// @copydoc SolverInterface::initialize
  void initialize();

  /// @copydoc SolverInterface::advance
  void advance(double computedTimeStepSize);

  /// @copydoc SolverInterface::finalize
  void finalize();

  ///@}

  ///@name Status Queries
  ///@{

  /// @copydoc SolverInterface::getMeshDimensions
  int getMeshDimensions(std::string_view meshName) const;

  /// @copydoc SolverInterface::getDataDimensions
  int getDataDimensions(std::string_view meshName, std::string_view dataName) const;

  /// @copydoc SolverInterface::isCouplingOngoing
  bool isCouplingOngoing() const;

  /// @copydoc SolverInterface::isTimeWindowComplete
  bool isTimeWindowComplete() const;

  /// @copydoc SolverInterface::getMaxTimeStepSize
  double getMaxTimeStepSize() const;

  ///@}

  ///@name Requirements
  ///@{

  /// @copydoc SolverInterface::requiresInitialData
  bool requiresInitialData();

  /// @copydoc SolverInterface::requiresReadingCheckpoint
  bool requiresReadingCheckpoint();

  /// @copydoc SolverInterface::requiresWritingCheckpoint
  bool requiresWritingCheckpoint();

  ///@}

  ///@name Mesh Access
  ///@anchor precice-mesh-access
  ///@{

  /// @copydoc SolverInterface::resetMesh
  void resetMesh(std::string_view meshName);

  /// @copydoc SolverInterface::hasMesh
  bool hasMesh(std::string_view meshName) const;

  /// @copydoc SolverInterface::requiresMeshConnectivityFor
  bool requiresMeshConnectivityFor(std::string_view meshName) const;

  /// @copydoc SolverInterface::requiresGradientDataFor
  bool requiresGradientDataFor(std::string_view meshName,
                               std::string_view dataName) const;

  /// @copydoc SolverInterface::setMeshVertex
  int setMeshVertex(
      std::string_view              meshName,
      ::precice::span<const double> position);

  /// @copydoc SolverInterface::getMeshVertexSize
  int getMeshVertexSize(std::string_view meshName) const;

  /// @copydoc SolverInterface::setMeshVertices
  void setMeshVertices(
      std::string_view              meshName,
      ::precice::span<const double> positions,
      ::precice::span<VertexID>     ids);

  /// @copydoc SolverInterface::setMeshEdge
  void setMeshEdge(
      std::string_view meshName,
      int              firstVertexID,
      int              secondVertexID);

  /// @copydoc SolverInterface::setMeshEdges
  void setMeshEdges(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  /// @copydoc SolverInterface::setMeshTriangle
  void setMeshTriangle(
      std::string_view meshName,
      int              firstVertexID,
      int              secondVertexID,
      int              thirdVertexID);

  /// @copydoc SolverInterface::setMeshTriangles
  void setMeshTriangles(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  /// @copydoc SolverInterface::setMeshQuad
  void setMeshQuad(
      std::string_view meshName,
      int              firstVertexID,
      int              secondVertexID,
      int              thirdVertexID,
      int              fourthVertexID);

  /// @copydoc SolverInterface::setMeshQuads
  void setMeshQuads(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  /// @copydoc SolverInterface::setMeshTetrahedron
  void setMeshTetrahedron(
      std::string_view meshName,
      int              firstVertexID,
      int              secondVertexID,
      int              thirdVertexID,
      int              fourthVertexID);

  /// @copydoc SolverInterface::setMeshTetrahedra
  void setMeshTetrahedra(
      std::string_view                meshName,
      ::precice::span<const VertexID> vertices);

  ///@}

  ///@name Data Access
  ///@{

  /// @copydoc SolverInterface::hasData
  bool hasData(
      std::string_view meshName,
      std::string_view dataName) const;

  /// @copydoc SolverInterface::readData
  void readData(
      std::string_view                meshName,
      std::string_view                dataName,
      ::precice::span<const VertexID> vertices,
      double                          relativeReadTime,
      ::precice::span<double>         values) const;

  /// @copydoc SolverInterface::writeData
  void writeData(
      std::string_view                meshName,
      std::string_view                dataName,
      ::precice::span<const VertexID> vertices,
      ::precice::span<const double>   values);

  /// @copydoc SolverInterface::writeGradientData
  void writeGradientData(
      std::string_view                meshName,
      std::string_view                dataName,
      ::precice::span<const VertexID> vertices,
      ::precice::span<const double>   gradients);

  ///@}

  /** @name Experimental Data Access
   * These API functions are \b experimental and may change in future versions.
   */
  ///@{

  /// @copydoc SolverInterface::setMeshAccessRegion
  void setMeshAccessRegion(std::string_view              meshName,
                           ::precice::span<const double> boundingBox) const;

  /// @copydoc SolverInterface::getMeshVerticesAndIDs
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
  SolverInterfaceImpl(SolverInterfaceImpl const &) = delete;

  /// Disable assignment construction
  SolverInterfaceImpl &operator=(SolverInterfaceImpl const &) = delete;

  /// Disable move construction
  SolverInterfaceImpl(SolverInterfaceImpl &&) = delete;

  /// Disable move assignment
  SolverInterfaceImpl &operator=(SolverInterfaceImpl &&) = delete;

private:
  mutable logging::Logger _log{"impl::SolverInterfaceImpl"};

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

  /// Represents the various states a SolverInterface can be in.
  enum struct State {
    Constructed, // Initial state of SolverInterface
    Initialized, // SolverInterface.initialize() triggers transition from State::Constructed to State::Initialized; mandatory
    Finalized    // SolverInterface.finalize() triggers transition form State::Initialized to State::Finalized; mandatory
  };

  /// Are experimental API calls allowed?
  bool _allowsExperimental = false;

  /// setMeshAccessRegion may only be called once
  mutable bool _accessRegionDefined = false;

  /// The current State of the solverinterface
  State _state{State::Constructed};

  /// Counts calls to advance for plotting.
  long int _numberAdvanceCalls = 0;

  /**
   * @brief Configures the coupling interface from the given xml file.
   *
   * Only after the configuration a reasonable state of a SolverInterfaceImpl
   * object is achieved.
   *
   * @param[in] configurationFileName Name (with path) of the xml config. file.
   */
  void configure(std::string_view configurationFileName);

  /**
   * @brief Configures the coupling interface with a prepared configuration.
   *
   * Can be used to configure the SolverInterfaceImpl without xml file. Requires
   * to manually setup the configuration object.
   */
  void configure(const config::SolverInterfaceConfiguration &configuration);

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

  /// Resets written data, displacements and mesh neighbors to export.
  void resetWrittenData();

  /// Determines participant accessing this interface from the configuration.
  impl::PtrParticipant determineAccessingParticipant(
      const config::SolverInterfaceConfiguration &config);

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

  /// To allow white box tests.
  friend struct Integration::Serial::Whitebox::TestConfigurationPeano;
  friend struct Integration::Serial::Whitebox::TestConfigurationComsol;
};

} // namespace impl
} // namespace precice
