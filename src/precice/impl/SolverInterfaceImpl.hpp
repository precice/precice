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

  /// @copydoc SolverInterface::isTimeWindowComplete
  bool isTimeWindowComplete() const;

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
  void resetMesh(MeshID meshID);

  /// @copydoc SolverInterface::hasMesh
  bool hasMesh(const std::string &meshName) const;

  /// @copydoc SolverInterface::hasMesh
  int getMeshID(const std::string &meshName) const;

  /// @copydoc SolverInterface::requiresMeshConnectivityFor
  bool requiresMeshConnectivityFor(int meshID) const;

  /// @copydoc SolverInterface::requiresGradientDataFor
  bool requiresGradientDataFor(int dataID) const;

  /// @copydoc SolverInterface::setMeshVertex
  int setMeshVertex(
      int           meshID,
      const double *position);

  /// @copydoc SolverInterface::getMeshVertexSize
  int getMeshVertexSize(MeshID meshID) const;

  /// @copydoc SolverInterface::setMeshVertices
  void setMeshVertices(
      int           meshID,
      int           size,
      const double *positions,
      int *         ids);

  /// @copydoc SolverInterface::setMeshEdge
  void setMeshEdge(
      MeshID meshID,
      int    firstVertexID,
      int    secondVertexID);

  /// @copydoc SolverInterface::setMeshEdges
  void setMeshEdges(
      int        meshID,
      int        size,
      const int *vertices);

  /// @copydoc SolverInterface::setMeshTriangle
  void setMeshTriangle(
      MeshID meshID,
      int    firstVertexID,
      int    secondVertexID,
      int    thirdVertexID);

  /// @copydoc SolverInterface::setMeshTriangles
  void setMeshTriangles(
      int        meshID,
      int        size,
      const int *vertices);

  /// @copydoc SolverInterface::setMeshQuad
  void setMeshQuad(
      MeshID meshID,
      int    firstVertexID,
      int    secondVertexID,
      int    thirdVertexID,
      int    fourthVertexID);

  /// @copydoc SolverInterface::setMeshQuads
  void setMeshQuads(
      int        meshID,
      int        size,
      const int *vertices);

  /// @copydoc SolverInterface::setMeshTetrahedron
  void setMeshTetrahedron(
      MeshID meshID,
      int    firstVertexID,
      int    secondVertexID,
      int    thirdVertexID,
      int    fourthVertexID);

  /// @copydoc SolverInterface::setMeshTetrahedra
  void setMeshTetrahedra(
      int        meshID,
      int        size,
      const int *vertices);

  ///@}

  ///@name Data Access
  ///@{

  /// @copydoc SolverInterface::hasData
  bool hasData(const std::string &dataName, MeshID meshID) const;

  /// @copydoc SolverInterface::getDataID
  int getDataID(const std::string &dataName, MeshID meshID) const;

  /// @copydoc SolverInterface::writeBlockVectorData
  void writeBlockVectorData(
      int           fromDataID,
      int           size,
      const int *   valueIndices,
      const double *values);

  /// @copydoc precice::SolverInterface::writeBlockVectorGradientData
  void writeBlockVectorGradientData(
      int           fromDataID,
      int           size,
      const int *   valueIndices,
      const double *gradientValues);

  /// @copydoc SolverInterface::writeVectorData
  void writeVectorData(
      int           fromDataID,
      int           valueIndex,
      const double *value);

  /// @copydoc precice::SolverInterface::writeVectorGradientData
  void writeVectorGradientData(
      int           fromDataID,
      int           valueIndex,
      const double *gradientValues);

  /// @copydoc SolverInterface::writeBlockScalarData
  void writeBlockScalarData(
      int           fromDataID,
      int           size,
      const int *   valueIndices,
      const double *values);

  /// @copydoc precice::SolverInterface::writeBlockScalarGradientData
  void writeBlockScalarGradientData(
      int           fromDataID,
      int           size,
      const int *   valueIndices,
      const double *gradientValues);

  /// @copydoc SolverInterface::writeScalarData
  void writeScalarData(
      int    fromDataID,
      int    valueIndex,
      double value);

  /// @copydoc precice::SolverInterface::writeScalarGradientData
  void writeScalarGradientData(
      int           fromDataID,
      int           valueIndex,
      const double *gradientValues);

  /// @copydoc SolverInterface::readBlockVectorData(int, int, const int*, double*) const
  void readBlockVectorData(
      int        toDataID,
      int        size,
      const int *valueIndices,
      double *   values) const;

  /// @copydoc SolverInterface::readBlockVectorData(int, int, const int*, double, double*) const
  void readBlockVectorData(
      int        toDataID,
      int        size,
      const int *valueIndices,
      double     relativeReadTime,
      double *   values) const;

  /// @copydoc SolverInterface::readVectorData(int, int, double*) const
  void readVectorData(
      int     toDataID,
      int     valueIndex,
      double *value) const;

  /// @copydoc SolverInterface::readVectorData(int, int, double, double*) const
  void readVectorData(
      int     toDataID,
      int     valueIndex,
      double  relativeReadTime,
      double *value) const;

  /// @copydoc SolverInterface::readBlockScalarData(int, int, const int*, double*) const
  void readBlockScalarData(
      int        toDataID,
      int        size,
      const int *valueIndices,
      double *   values) const;

  /// @copydoc SolverInterface::readBlockScalarData(int, int, const int*, double, double*) const
  void readBlockScalarData(
      int        toDataID,
      int        size,
      const int *valueIndices,
      double     relativeReadTime,
      double *   values) const;

  /// @copydoc SolverInterface::readScalarData(int, int, double&) const
  void readScalarData(
      int     toDataID,
      int     valueIndex,
      double &value) const;

  /// @copydoc SolverInterface::readScalarData(int, int, double, double&) const
  void readScalarData(
      int     toDataID,
      int     valueIndex,
      double  relativeReadTime,
      double &value) const;

  ///@}

  /** @name Experimental Data Access
   * These API functions are \b experimental and may change in future versions.
   */
  ///@{

  /// @copydoc SolverInterface::setMeshAccessRegion
  void setMeshAccessRegion(const int     meshID,
                           const double *boundingBox) const;

  /// @copydoc SolverInterface::getMeshVerticesAndIDs
  void getMeshVerticesAndIDs(
      const int meshID,
      const int size,
      int *     ids,
      double *  coordinates) const;

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
  /**
   * @brief Generic constructor for SolverInterfaceImpl.
   *
   * Use the parameter communicator to specify a custom global MPI communicator.
   * Pass a null pointer to signal preCICE to use MPI_COMM_WORLD.
   *
   * @param[in] participantName Name of the participant using the interface. Has to
   *        match the name given for a participant in the xml configuration file.
   * @param[in] configurationFileName Name (with path) of the xml configuration file.
   * @param[in] solverProcessIndex If the solver code runs with several processes,
   *        each process using preCICE has to specify its index, which has to start
   *        from 0 and end with solverProcessSize - 1.
   * @param[in] solverProcessSize The number of solver processes using preCICE.
   * @param[in] communicator A pointer to an MPI_Comm to use as communicator.
   * @param[in] allowNullptr    Accept nullptr for communicator.
   */
  SolverInterfaceImpl(
      std::string        participantName,
      const std::string &configurationFileName,
      int                solverProcessIndex,
      int                solverProcessSize,
      void *             communicator,
      bool               allowNullptr);

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
      int        toDataID,
      int        size,
      const int *valueIndices,
      double     relativeReadTime,
      double *   values) const;

  void readVectorDataImpl(
      int     toDataID,
      int     valueIndex,
      double  relativeReadTime,
      double *value) const;

  void readBlockScalarDataImpl(
      int        toDataID,
      int        size,
      const int *valueIndices,
      double     relativeReadTime,
      double *   values) const;

  void readScalarDataImpl(
      int     toDataID,
      int     valueIndex,
      double  relativeReadTime,
      double &value) const;

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
