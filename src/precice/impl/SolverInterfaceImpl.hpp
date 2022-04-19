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

/// Implementation of solver interface.
class SolverInterfaceImpl {
public:
  /**
   * @copydoc SolverInterface::SolverInterface(std::string, const std::string, int, int)
   */
  SolverInterfaceImpl(
      std::string        participantName,
      const std::string &configurationFileName,
      int                solverProcessIndex,
      int                solverProcessSize);

  /// Deleted copy constructor
  SolverInterfaceImpl(SolverInterfaceImpl const &) = delete;

  /// Deleted copy assignment
  SolverInterfaceImpl &operator=(SolverInterfaceImpl const &) = delete;

  /// Deleted move constructor
  SolverInterfaceImpl(SolverInterfaceImpl &&) = delete;

  /// Deleted move assignment
  SolverInterfaceImpl &operator=(SolverInterfaceImpl &&) = delete;

  /**
   * @copydoc SolverInterface::SolverInterface(std::string, const std::string, int, int, void*)
   */
  SolverInterfaceImpl(
      std::string        participantName,
      const std::string &configurationFileName,
      int                solverProcessIndex,
      int                solverProcessSize,
      void *             communicator);

  /** Ensures that finalize() has been called.
   *
   * @see finalize()
   */
  ~SolverInterfaceImpl();

  /// @copydoc SolverInterface::initialize()
  double initialize();

  /// @copydoc SolverInterface::initializeData()
  void initializeData();

  /// @copydoc SolverInterface::advance(double)
  double advance(double computedTimestepLength);

  /// @copydoc SolverInterface::finalize()
  void finalize();

  /// @copydoc SolverInterface::getDimensions()
  int getDimensions() const;

  /// @copydoc SolverInterface::isCouplingOngoing()
  bool isCouplingOngoing() const;

  /// @copydoc SolverInterface::isReadDataAvailable()
  bool isReadDataAvailable() const;

  /// @copydoc SolverInterface::isWriteDataRequired(double)
  bool isWriteDataRequired(double computedTimestepLength) const;

  /// @copydoc SolverInterface::isTimeWindowComplete()
  bool isTimeWindowComplete() const;

  /// @copydoc SolverInterface::hasToEvaluateSurrogateModel()
  bool hasToEvaluateSurrogateModel() const;

  /// @copydoc SolverInterface::hasToEvaluateFineModel()
  bool hasToEvaluateFineModel() const;

  /// @copydoc SolverInterface::isActionRequired(std::string)
  bool isActionRequired(const std::string &action) const;

  /// @copydoc SolverInterface::markActionFulfilled(std::string)
  void markActionFulfilled(const std::string &action);

  /// @copydoc SolverInterface::hasMesh(std::string)
  bool hasMesh(const std::string &meshName) const;

  /// @copydoc SolverInterface::resetMesh(int)
  void resetMesh(MeshID meshID);

  /// @copydoc SolverInterface::hasMesh(std::string)
  int getMeshID(const std::string &meshName) const;

  /// @copydoc SolverInterface::getMeshIDs()
  std::set<int> getMeshIDs() const;

  /// @copydoc SolverInterface::isMeshConnectivityRequired()
  bool isMeshConnectivityRequired(int meshID) const;

  /// @copydoc SolverInterface::setMeshVertex(int, double*)
  int setMeshVertex(
      int           meshID,
      const double *position);

  /// @copydoc SolverInterface::getMeshVertexSize(int)
  int getMeshVertexSize(MeshID meshID) const;

  /// @copydoc SolverInterface::setMeshVertices(int, int, double*, int*)
  void setMeshVertices(
      int           meshID,
      int           size,
      const double *positions,
      int *         ids);

  /// @copydoc SolverInterface::getMeshVertices(int, int, int*, double*)
  void getMeshVertices(
      int        meshID,
      size_t     size,
      const int *ids,
      double *   positions) const;

  /// @copydoc SolverInterface::getMeshVertexIDsFromPositions(int, int, double*, int*)
  void getMeshVertexIDsFromPositions(
      int           meshID,
      size_t        size,
      const double *positions,
      int *         ids) const;

  /// @copydoc SolverInterface::setMeshEdge(int, int, int)
  int setMeshEdge(
      MeshID meshID,
      int    firstVertexID,
      int    secondVertexID);

  /// @copydoc SolverInterface::setMeshTriangle(int, int, int, int)
  void setMeshTriangle(
      MeshID meshID,
      int    firstEdgeID,
      int    secondEdgeID,
      int    thirdEdgeID);

  /// @copydoc SolverInterface::setMeshTriangleWithEdges(int, int, int, int)
  void setMeshTriangleWithEdges(
      MeshID meshID,
      int    firstVertexID,
      int    secondVertexID,
      int    thirdVertexID);

  /// @copydoc SolverInterface::setMeshQuad(int, int, int, int, int)
  void setMeshQuad(
      MeshID meshID,
      int    firstEdgeID,
      int    secondEdgeID,
      int    thirdEdgeID,
      int    fourthEdgeID);

  /// @copydoc SolverInterface::setMeshQuadWithEdges(int, int, int, int, int)
  void setMeshQuadWithEdges(
      MeshID meshID,
      int    firstVertexID,
      int    secondVertexID,
      int    thirdVertexID,
      int    fourthVertexID);

  /// @copydoc SolverInterface::hasData(std::string, int)
  bool hasData(const std::string &dataName, MeshID meshID) const;

  /// @copydoc SolverInterface::getDataID(std::string, int)
  int getDataID(const std::string &dataName, MeshID meshID) const;

  /// @copydoc SolverInterface::mapWriteDataFrom(int)
  void mapWriteDataFrom(int fromMeshID);

  /// @copydoc SolverInterface::mapReadDataTo(int)
  void mapReadDataTo(int toMeshID);

  /// @copydoc SolverInterface::writeBlockVectorData(int, int, const int*, const double*)
  void writeBlockVectorData(
      int           fromDataID,
      int           size,
      const int *   valueIndices,
      const double *values);

  /// @copydoc SolverInterface::writeVectorData(int, int, const double*)
  void writeVectorData(
      int           fromDataID,
      int           valueIndex,
      const double *value);

  /// @copydoc SolverInterface::writeBlockScalarData(int, int, const int*, const double*)
  void writeBlockScalarData(
      int           fromDataID,
      int           size,
      const int *   valueIndices,
      const double *values);

  /// @copydoc SolverInterface::writeScalarData(int, int, double)
  void writeScalarData(
      int    fromDataID,
      int    valueIndex,
      double value);

  /// @copydoc SolverInterface::readBlockVectorData(int, int, const int*, double*)
  void readBlockVectorData(
      int        toDataID,
      int        size,
      const int *valueIndices,
      double *   values) const;

  /// @copydoc SolverInterface::readVectorData(int, int, double*)
  void readVectorData(
      int     toDataID,
      int     valueIndex,
      double *value) const;

  /// @copydoc SolverInterface::readBlockScalarData(int, int, const int*, double*)
  void readBlockScalarData(
      int        toDataID,
      int        size,
      const int *valueIndices,
      double *   values) const;

  /// @copydoc SolverInterface::readScalarData(int, int, double&)
  void readScalarData(
      int     toDataID,
      int     valueIndex,
      double &value) const;

  /// @copydoc SolverInterface::setMeshAccessRegion(const int, const double*)
  void setMeshAccessRegion(const int     meshID,
                           const double *boundingBox) const;

  /// @copydoc SolverInterface::getMeshVerticesAndIDs(const int, const int, int*, double*)
  void getMeshVerticesAndIDs(
      const int meshID,
      const int size,
      int *     ids,
      double *  coordinates) const;

  /// @copydoc SolverInterface::exportMesh(const std::string)
  void exportMesh(const std::string &filenameSuffix) const;

  /// Allows to access a registered mesh
  /// @todo make private and use fixture.
  const mesh::Mesh &mesh(const std::string &meshName) const;

private:
  mutable logging::Logger _log{"impl::SolverInterfaceImpl"};

  std::string _accessorName;

  int _accessorProcessRank;

  int _accessorCommunicatorSize;

  impl::PtrParticipant _accessor;

  /// Spatial dimensions of problem.
  int _dimensions = 0;

  utils::MultiLock<int> _meshLock;

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

  /// Initializes communication between master and slaves.
  void initializeMasterSlaveCommunication();

  /// Syncs the timestep between slaves and master (all timesteps should be the same!)
  void syncTimestep(double computedTimestepLength);

  /// Which channels to close in closeCommunicationChannels()
  enum class CloseChannels : bool {
    All         = false,
    Distributed = true
  };

  /// Syncs the masters of all connected participants
  void closeCommunicationChannels(CloseChannels cc);

  /// To allow white box tests.
  friend struct Integration::Serial::Whitebox::TestConfigurationPeano;
  friend struct Integration::Serial::Whitebox::TestConfigurationComsol;
};

} // namespace impl
} // namespace precice
