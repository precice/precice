#pragma once

#include <map>
#include <set>
#include <stddef.h>
#include <string>
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
   * @copydoc SolverInterface::SolverInterface(const std::string&, const std::string&, int, int)
   */
  SolverInterfaceImpl(
      std::string        participantName,
      const std::string &configurationFileName,
      int                solverProcessIndex,
      int                solverProcessSize);

  /**
   * @copybrief SolverInterface::SolverInterface(const std::string&, const std::string&, int, int, void*)
   *
   * Use the parameter communicator to specify a custom global MPI communicator.
   * Pass a null pointer to signal preCICE to use MPI_COMM_WORLD.
   *
   * @copydetails SolverInterface::SolverInterface(const std::string&, const std::string&, int, int, void*)
   */
  SolverInterfaceImpl(
      std::string        participantName,
      const std::string &configurationFileName,
      int                solverProcessIndex,
      int                solverProcessSize,
      void *             communicator);

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
  double initialize();

  /// @copydoc SolverInterface::initializeData
  void initializeData();

  /// @copydoc SolverInterface::advance
  double advance(double computedTimestepLength);

  /// @copydoc SolverInterface::finalize
  void finalize();

  ///@}

  ///@name Status Queries
  ///@{

  /// @copydoc SolverInterface::getDimensions
  int getDimensions() const;

  /// @copydoc SolverInterface::isCouplingOngoing
  bool isCouplingOngoing() const;

  /// @copydoc SolverInterface::isReadDataAvailable
  bool isReadDataAvailable() const;

  /// @copydoc SolverInterface::isWriteDataRequired
  bool isWriteDataRequired(double computedTimestepLength) const;

  /// @copydoc SolverInterface::isTimeWindowComplete
  bool isTimeWindowComplete() const;

  /// @copydoc SolverInterface::hasToEvaluateSurrogateModel
  bool hasToEvaluateSurrogateModel() const;

  /// @copydoc SolverInterface::hasToEvaluateFineModel
  bool hasToEvaluateFineModel() const;

  ///@}

  ///@name Action Methods
  ///@{

  /// @copydoc SolverInterface::isActionRequired
  bool isActionRequired(const std::string &action) const;

  /// @copydoc SolverInterface::markActionFulfilled
  void markActionFulfilled(const std::string &action);

  ///@}

  ///@name Mesh Access
  ///@anchor precice-mesh-access
  ///@{

  /// @copydoc SolverInterface::resetMesh
  void resetMesh(MeshID meshID);

  /// @copydoc SolverInterface::hasMesh
  bool hasMesh(const std::string &meshName) const;

  /// @copydoc SolverInterface::hasMesh
  MeshID getMeshID(const std::string &meshName) const;

  /// @copydoc SolverInterface::getMeshIDs
  std::set<MeshID> getMeshIDs() const;

  /// @copydoc SolverInterface::isMeshConnectivityRequired
  bool isMeshConnectivityRequired(MeshID meshID) const;

  /// @copydoc SolverInterface::isGradientDataRequired
  bool isGradientDataRequired(DataID dataID) const;

  /// @copydoc SolverInterface::setMeshVertex
  VertexID setMeshVertex(
      MeshID        meshID,
      const double *position);

  /// @copydoc SolverInterface::getMeshVertexSize
  Size getMeshVertexSize(MeshID meshID) const;

  /// @copydoc SolverInterface::setMeshVertices
  void setMeshVertices(
      MeshID        meshID,
      Size          size,
      const double *positions,
      VertexID *    ids);

  /// @copydoc SolverInterface::getMeshVertices
  void getMeshVertices(
      MeshID          meshID,
      Size            size,
      const VertexID *ids,
      double *        positions) const;

  /// @copydoc SolverInterface::getMeshVertexIDsFromPositions
  void getMeshVertexIDsFromPositions(
      MeshID        meshID,
      Size          size,
      const double *positions,
      VertexID *    ids) const;

  /// @copydoc SolverInterface::setMeshEdge
  EdgeID setMeshEdge(
      MeshID   meshID,
      VertexID firstVertexID,
      VertexID secondVertexID);

  /// @copydoc SolverInterface::setMeshTriangle
  void setMeshTriangle(
      MeshID meshID,
      EdgeID firstEdgeID,
      EdgeID secondEdgeID,
      EdgeID thirdEdgeID);

  /// @copydoc SolverEdgeIDerface::setMeshTriangleWithEdges
  void setMeshTriangleWithEdges(
      MeshID   meshID,
      VertexID firstVertexID,
      VertexID secondVertexID,
      VertexID thirdVertexID);

  /// @copydoc SolverEdgeIDerface::setMeshQuad
  void setMeshQuad(
      MeshID meshID,
      EdgeID firstEdgeID,
      EdgeID secondEdgeID,
      EdgeID thirdEdgeID,
      EdgeID fourthEdgeID);

  /// @copydoc SolverInterface::setMeshQuadWithEdges
  void setMeshQuadWithEdges(
      MeshID   meshID,
      VertexID firstVertexID,
      VertexID secondVertexID,
      VertexID thirdVertexID,
      VertexID fourthVertexID);

  ///@}

  ///@name Data Access
  ///@{

  /// @copydoc SolverInterface::hasData
  bool hasData(const std::string &dataName, MeshID meshID) const;

  /// @copydoc SolverInterface::getDataID
  DataID getDataID(const std::string &dataName, MeshID meshID) const;

  /// @copydoc SolverInterface::mapWriteDataFrom
  void mapWriteDataFrom(MeshID fromMeshID);

  /// @copydoc SolverInterface::mapReadDataTo
  void mapReadDataTo(MeshID toMeshID);

  /// @copydoc SolverInterface::writeBlockVectorData
  void writeBlockVectorData(
      DataID          fromDataID,
      Size            size,
      const VertexID *valueIndices,
      const double *  values);

  /// @copydoc precice::SolverInterface::writeBlockVectorGradientData
  void writeBlockVectorGradientData(
      DataID          fromDataID,
      Size            size,
      const VertexID *valueIndices,
      const double *  gradientValues,
      bool            rowsFirst = false);

  /// @copydoc SolverInterface::writeVectorData
  void writeVectorData(
      DataID        fromDataID,
      VertexID      valueIndex,
      const double *value);

  /// @copydoc precice::SolverInterface::writeVectorGradientData
  void writeVectorGradientData(
      DataID        fromDataID,
      VertexID      valueIndex,
      const double *gradientValues,
      bool          rowsFirst = false);

  /// @copydoc SolverInterface::writeBlockScalarData
  void writeBlockScalarData(
      DataID          fromDataID,
      Size            size,
      const VertexID *valueIndices,
      const double *  values);

  /// @copydoc precice::SolverInterface::writeBlockScalarGradientData
  void writeBlockScalarGradientData(
      DataID          fromDataID,
      Size            size,
      const VertexID *valueIndices,
      const double *  gradientValues);

  /// @copydoc SolverInterface::writeScalarData
  void writeScalarData(
      DataID   fromDataID,
      VertexID valueIndex,
      double   value);

  /// @copydoc precice::SolverInterface::writeScalarGradientData
  void writeScalarGradientData(
      DataID        fromDataID,
      VertexID      valueIndex,
      const double *gradientValues);

  /// @copydoc SolverInterface::readBlockVectorData(int, int, const int*, double*) const
  void readBlockVectorData(
      DataID          toDataID,
      Size            size,
      const VertexID *valueIndices,
      double *        values) const;

  /// @copydoc SolverInterface::readBlockVectorData(int, int, const int*, double, double*) const
  void readBlockVectorData(
      DataID          toDataID,
      Size            size,
      const VertexID *valueIndices,
      double          relativeReadTime,
      double *        values) const;

  /// @copydoc SolverInterface::readVectorData(int, int, double*) const
  void readVectorData(
      DataID   toDataID,
      VertexID valueIndex,
      double * value) const;

  /// @copydoc SolverInterface::readVectorData(int, int, double, double*) const
  void readVectorData(
      DataID   toDataID,
      VertexID valueIndex,
      double   relativeReadTime,
      double * value) const;

  /// @copydoc SolverInterface::readBlockScalarData(int, int, const int*, double*) const
  void readBlockScalarData(
      DataID          toDataID,
      Size            size,
      const VertexID *valueIndices,
      double *        values) const;

  /// @copydoc SolverInterface::readBlockScalarData(int, int, const int*, double, double*) const
  void readBlockScalarData(
      DataID          toDataID,
      Size            size,
      const VertexID *valueIndices,
      double          relativeReadTime,
      double *        values) const;

  /// @copydoc SolverInterface::readScalarData(int, int, double&) const
  void readScalarData(
      DataID   toDataID,
      VertexID valueIndex,
      double & value) const;

  /// @copydoc SolverInterface::readScalarData(int, int, double, double&) const
  void readScalarData(
      DataID   toDataID,
      VertexID valueIndex,
      double   relativeReadTime,
      double & value) const;

  ///@}

  /** @name Experimental Data Access
   * These API functions are \b experimental and may change in future versions.
   */
  ///@{

  /// @copydoc SolverInterface::setMeshAccessRegion
  void setMeshAccessRegion(const MeshID  meshID,
                           const double *boundingBox) const;

  /// @copydoc SolverInterface::getMeshVerticesAndIDs
  void getMeshVerticesAndIDs(
      const MeshID meshID,
      const Size   size,
      VertexID *   ids,
      double *     coordinates) const;

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

  utils::MultiLock<MeshID> _meshLock;

  /// mesh name to mesh ID mapping.
  std::map<std::string, MeshID> _meshIDs;

  std::map<std::string, m2n::BoundM2N> _m2ns;

  /// Holds information about solvers participating in the coupled simulation.
  std::vector<impl::PtrParticipant> _participants;

  cplscheme::PtrCouplingScheme _couplingScheme;

  /// Represents the various states a SolverInterface can be in.
  enum struct State {
    Constructed, // Initial state of SolverInterface
    Initialized, // SolverInterface.initialize() triggers transition from State::Constructed to State::Initialized; mandatory
    Finalized    // SolverInterface.finalize() triggers transition form State::Initialized or State::InitializedData to State::Finalized; mandatory
  };

  /// SolverInterface.initializeData() triggers transition from false to true.
  bool _hasInitializedData = false;

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
  void configure(const std::string &configurationFileName);

  /**
   * @brief Configures the coupling interface with a prepared configuration.
   *
   * Can be used to configure the SolverInterfaceImpl without xml file. Requires
   * to manually setup the configuration object.
   */
  void configure(const config::SolverInterfaceConfiguration &configuration);

  void configureM2Ns(const m2n::M2NConfiguration::SharedPointer &config);

  /// Implementation of read functions.
  void readBlockVectorDataImpl(
      DataID          toDataID,
      Size            size,
      const VertexID *valueIndices,
      double          relativeReadTime,
      double *        values) const;

  void readVectorDataImpl(
      DataID   toDataID,
      VertexID valueIndex,
      double   relativeReadTime,
      double * value) const;

  void readBlockScalarDataImpl(
      DataID          toDataID,
      Size            size,
      const VertexID *valueIndices,
      double          relativeReadTime,
      double *        values) const;

  void readScalarDataImpl(
      DataID   toDataID,
      VertexID valueIndex,
      double   relativeReadTime,
      double & value) const;

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
  void computeMappings(const utils::ptr_vector<MappingContext> &contexts, const std::string &mappingType);

  /// Helper for mapWrittenData and mapReadData
  void clearMappings(utils::ptr_vector<MappingContext> contexts);

  /// Computes, performs, and resets all suitable write mappings.
  void mapWrittenData();

  /// Computes, performs, and resets all suitable read mappings.
  void mapReadData();

  /**
   * @brief Performs all data actions with given timing.
   *
   * @param[in] timings the timings of the action.
   * @param[in] time the current total simulation time.
   * @param[in] timeStepSize Length of last time step computed.
   * @param[in] computedTimeWindowPart Sum of all time steps within current time window, i.e. part that is already computed.
   * @param[in] timeWindowSize Current time window size.
   */
  void performDataActions(
      const std::set<action::Action::Timing> &timings,
      double                                  time,
      double                                  timeStepSize,
      double                                  computedTimeWindowPart,
      double                                  timeWindowSize);

  /// Resets written data, displacements and mesh neighbors to export.
  void resetWrittenData();

  /// Determines participant accessing this interface from the configuration.
  impl::PtrParticipant determineAccessingParticipant(
      const config::SolverInterfaceConfiguration &config);

  /// Initializes intra-participant communication.
  void initializeIntraCommunication();

  /// Syncs the timestep between all ranks (all timesteps should be the same!)
  void syncTimestep(double computedTimestepLength);

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
