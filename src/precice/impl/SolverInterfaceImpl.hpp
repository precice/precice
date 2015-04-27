// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_SOLVERINTERFACEIMPL_HPP_
#define PRECICE_IMPL_SOLVERINTERFACEIMPL_HPP_

#include "precice/MeshHandle.hpp"
#include "precice/Constants.hpp"
#include "precice/ClosestMesh.hpp"
#include "precice/VoxelPosition.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "precice/impl/DataContext.hpp"
#include "action/Action.hpp"
#include "utils/Dimensions.hpp"
#include "boost/noncopyable.hpp"
#include "io/Constants.hpp"
#include "query/ExportVTKNeighbors.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "tarch/la/WrappedVector.h"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include <string>
#include <vector>
#include <set>

namespace precice {
  namespace impl {
    class RequestManager;
  }
  namespace config {
    class SolverInterfaceConfiguration;
  }
  namespace tests {
    class SolverInterfaceTest;
    class SolverInterfaceTestRemote;
    class SolverInterfaceTestGeometry;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace impl {

/**
 * @brief Implementation of solver interface.
 */
class SolverInterfaceImpl : private boost::noncopyable
{
public:

  /**
   * @brief Constructor.
   *
   * A solver that wants to use the SolverInterfaceImpl must instatiate an object
   * of this class. The object has to be configured by one of the configure
   * methods before it has a reasonable state and can be used.
   *
   * @param accessorName [IN] Name of the solver using the interface. Has to
   *                          match the name given for a participant in the
   *                          xml configuration file.
   */
  SolverInterfaceImpl (
    const std::string& accessorName,
    int                accessorProcessRank,
    int                accessorCommunicatorSize,
    bool               serverMode );

  /**
   * @brief Destructor.
   */
  ~SolverInterfaceImpl();

  /**
   * @brief Configures the coupling interface from the given xml file.
   *
   * Only after the configuration a reasonable state of a SolverInterfaceImpl
   * object is achieved.
   *
   * @param configurationFileName [IN] Name (with path) of the xml config. file.
   */
  void configure ( const std::string& configurationFileName );

  /**
   * @brief Configures the coupling interface with a prepared configuration.
   *
   * Can be used to configure the SolverInterfaceImpl without xml file. Requires
   * to manually setup the configuration object.
   */
  void configure ( const config::SolverInterfaceConfiguration& configuration );

  /**
   * @brief Initializes all coupling data and starts (coupled) simulation.
   *
   * - Initiates MPI communication, if not done yet.
   * - Sets up a connection to the other participants of the coupled simulation.
   * - Creates all meshes, solver meshes need to be submitted before.
   * - Receives first coupling data, when the solver is not starting the
   *   coupled simulation.
   * - Determines length of the first timestep to be computed.
   *
   * @return Maximum length of first timestep to be computed by the solver.
   */
  double initialize();

  /**
   * @brief Sets the sofar written data as initial value for the coupling scheme.
   *
   * Erases the written data afterwards.
   */
  void initializeData();

  /**
   * @brief Exchanges coupling data and advances coupling state.
   *
   * - Updates the topology of geometries changed by the solver.
   * - Sends and resets coupling data written by solver to coupling partners.
   * - Sends geometries with changed topology to coupling partners.
   * - Receives coupling data read by solver.
   * - Receives geometries that are changed by other solvers.
   * - Exchanges and computes information regarding the state of the coupled
   *   simulation.
   *
   * @param computedTimestepLength [IN] Length of timestep computed by solver.
   * @return Maximum length of next timestep to be computed by solver.
   */
  double advance ( double computedTimestepLength );

  /**
   * @brief Finalizes the coupled simulation.
   *
   * - Tears down communication means used for coupling.
   * - Finalizes MPI if initialized in initializeCoupling.
   */
  void finalize();

  /**
   * @brief Returns the number of spatial dimensions for the coupling.
   *
   * The number of dimensions is fixed on configuration.
   */
  int getDimensions() const;

  /**
   * @brief Returns true, if the coupled simulation is still ongoing.
   *
   * The information to decide about the continuation of the coupled simulation
   * is retreived in the function initializeCoupling and updated in the
   * function exchangeData.
   */
  bool isCouplingOngoing();

  /**
   * @brief Returns true, if new data to be read is available.
   */
  bool isReadDataAvailable();

  /**
   * @brief Returns true, if new data has to be written.
   */
  bool isWriteDataRequired ( double computedTimestepLength );

  /**
   * @brief Returns true, if a global timestep is completed.
   */
  bool isTimestepComplete();

  /**
   * @brief Returns true, if provided name of action is required.
   *
   * Some features of preCICE require a solver to perform specific actions, in
   * order to be in valid state for a coupled simulation. A solver is made
   * eligible to use those features, by querying for the required actions,
   * performing them on demand, and calling fulfilledAction() to signalize
   * preCICE the correct behavior of the solver.
   */
  bool isActionRequired (	const std::string& action );

  /**
   * @brief Tells preCICE that a required action has been fulfilled by a solver.
   *
   * For more details see method requireAction().
   */
  void fulfilledAction ( const std::string& action );

  /**
   * @brief Returns true, if the mesh with given name is used.
   */
  bool hasMesh ( const std::string& meshName ) const;

  /**
   * @brief Returns the geometry ID belonging to the geometry with given name.
   *
   * The existing names are determined from the configuration.
   */
  int getMeshID (	const std::string& meshName );

  /**
   * @brief Returns all mesh IDs (besides sub-ids).
   */
  std::set<int> getMeshIDs();

  /**
   * @brief Returns true, if the data with given name is used in the given mesh.
   */
  bool hasData ( const std::string& dataName, int meshID );

  /**
   * @brief Returns data id corresponding to the given name (from configuration)
   * and mesh.
   */
  int getDataID ( const std::string& dataName, int meshID );

  /**
   * @brief Find out position of point relative to geometries.
   *
   * This method utililizes the caching mechanism of spacetees, and runs much
   * faster for non-local queries than inquireClosestMesh().
   *
   * @param point [IN] Point for which the inquiry is done.
   * @param meshIDs [IN] Mesh IDs of the geometries to be considered.
   * @return Position ID (see precice::constants::positionXYZ()).
   */
  int inquirePosition (
    const double*        point,
    const std::set<int>& meshIDs );

  /**
   * @brief Find out, which geometry is closest to the specified point
   *
   * This function works only for a continuous geometry, not for data sets.
   * The distance is positive, if the inquired point lies out of the
   * geometry, and negative, if it lies within the geometry.
   *
   * @param inquiredPoint       [IN] Coordinates of the point inquired
   * @param distanceToGeometry [OUT] Distance to the next geometry
   * @param geometryID         [OUT] ID of the next geometry
   *
   * @return AdjacentGeometry object containing geometry id and distance
   */
  ClosestMesh inquireClosestMesh (
    const double*        point,
    const std::set<int>& meshIDs );

  /**
   * @brief Inquires, whether voxel is inside, outside or on geometry surface
   *
   * Objects touching the voxel surface are not considered to be in the voxel.
   *
   * @param inquiryCenter   [IN] Center point of voxel
   * @param halfLengthVoxel [IN] Half length of voxel sides
   *
   * @return Constant, specifying if voxel is inside, outside or on geometry
   */
  VoxelPosition inquireVoxelPosition (
    const double*        voxelCenter,
    const double*        voxelHalflengths,
    bool                 includeBoundaries,
//    bool                 storeContent,
    const std::set<int>& meshIDs );

  /**
   * @brief Resets mesh with given ID.
   *
   * Has to be called, everytime the positions for data to be mapped
   * changes. Only has an effect, if the mapping used is non-stationary and
   * non-incremental.
   */
  void resetMesh ( int meshID );

  /**
   * @brief Sets several spatial positions for a mesh.
   *
   * @param ids [OUT] IDs for data from given positions.
   */
  void setMeshVertices (
    int     meshID,
    int     size,
    double* positions,
    int*    ids );

  /**
   * @brief Gets spatial positions of vertices for given IDs.
   *
   * @param ids [IN] IDs obtained when setting write positions.
   * @param positions [OUT] Positions corresponding to IDs.
   */
  void getMeshVertices (
    int     meshID,
    size_t  size,
    int*    ids,
    double* positions );

  /**
   * @brief Gets vertex data ids from positions.
   *
   * @param size [IN] Number of positions, ids.
   * @param positions [IN] Positions (x,y,z,x,y,z,...) to find ids for.
   * @param ids [OUT] IDs corresponding to positions.
   */
  void getMeshVertexIDsFromPositions (
    int     meshID,
    size_t  size,
    double* positions,
    int*    ids );

  /**
   * @brief Returns the number of nodes of a mesh.
   */
  int getMeshVertexSize ( int meshID );

  /**
   * @brief Set the position of a solver mesh vertex.
   *
   * @return Vertex ID to be used when setting an edge.
   */
  int setMeshVertex (
    int           meshID,
    const double* position );

  /**
   * @brief Set an edge of a solver mesh.
   *
   * @return Index of the edge to be used when setting a triangle.
   */
  int setMeshEdge (
    int meshID,
    int firstVertexID,
    int secondVertexID );

  /**
   * @brief Set a triangle of a solver mesh.
   */
  void setMeshTriangle (
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID );

  /**
   * @brief Sets a triangle and creates/sets edges automatically of a solver mesh.
   */
  void setMeshTriangleWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID );

  /**
   * @brief Set a quadrangle of a solver mesh.
   */
  void setMeshQuad (
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID,
    int fourthEdgeID );

  /**
   * @brief Sets a quadrangle and creates/sets edges automatically of a solver mesh.
   */
  void setMeshQuadWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID );

  /**
   * @brief Computes and maps all write data mapped from mesh with given ID.
   *
   * Is automatically called in advance, if not called manually before.
   */
  void mapWriteDataFrom(int fromMeshID);

  /**
   * @brief Computes and maps all read data mapped to mesh with given ID.
   */
  void mapReadDataTo(int toMeshID);

  /**
   * @brief Writes vector data values given as block.
   *
   * The block must contain the vector values in the following form:
   * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
   * the number of vector values. In 2D, the z-components are removed.
   *
   * @param fromDataID [IN] ID of the data to be written.
   * @param size [IN] Number of valueIndices, and number of values * dimensions.
   * @param values [IN] Values of the data to be written.
   */
  void writeBlockVectorData (
    int     fromDataID,
    int     size,
    int*    valueIndices,
    double* values );


  /**
   * @brief Write vectorial data to the geometry interface
   *
   * The exact mapping and communication must be specified in XYZ.
   *
   * @param fromDataID       [IN] ID of the data to be written, e.g. 1 = forces
   * @param dataPosition [IN] Position (coordinate, e.g.) of data to be written
   * @param dataValue    [IN] Value of the data to be written
   */
  void writeVectorData (
    int           fromDataID,
    int           valueIndex,
    const double* value );

  /**
   * @brief Writes scalar data values given as block.
   *
   * @param fromDataID [IN] ID of the data to be written.
   * @param size [IN] Number of valueIndices, and number of values.
   * @param values [IN] Values of the data to be written.
   */
  void writeBlockScalarData (
    int     fromDataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Write scalar data to the geometry interface
   *
   * The exact mapping and communication must be specified in XYZ.
   *
   * @param fromDataID       [IN] ID of the data to be written (2 = temperature, e.g.)
   * @param dataPosition [IN] Position (coordinate, e.g.) of data to be written
   * @param dataValue    [IN] Value of the data to be written
   */
  void writeScalarData(
    int    fromDataID,
    int    valueIndex,
    double value );

  /**
   * @brief Reads vector data values given as block.
   *
   * The block contains the vector values in the following form:
   * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
   * the number of vector values. In 2D, the z-components are removed.
   *
   * @param toDataID [IN] ID of the data to be read.
   * @param size [IN] Number of indices, and number of values * dimensions.
   * @param valueIndices [IN] Indices (from setReadPosition()) of data values.
   * @param values [IN] Values of the data to be read.
   */
  void readBlockVectorData (
    int     toDataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Reads vector data from the coupling mesh.
   *
   * @param toDataID       [IN]  ID of the data to be read, e.g. 1 = forces
   * @param dataPosition [IN]  Position (coordinate, e.g.) of data to be read
   * @param dataValue    [OUT] Read data value
   */
  void readVectorData (
    int     toDataID,
    int     valueIndex,
    double* value );

  /**
   * @brief Reads scalar data values given as block.
   *
   * @param toDataID [IN] ID of the data to be written.
   * @param size [IN] Number of valueIndices, and number of values.
   * @param values [IN] Values of the data to be written.
   */
  void readBlockScalarData (
    int     toDataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Read scalar data from the geometry interface.
   *
   * The exact mapping and communication must be specified in XYZ.
   *
   * @param toDataID       [IN]  ID of the data to be read, e.g. 2 = temperatures
   * @param dataPosition [IN]  Position (coordinate, e.g.) of data to be read
   * @param dataValue    [OUT] Read data value
   */
  void readScalarData (
    int     toDataID,
    int     valueIndex,
    double& value );

  /**
   * @brief Sets the location for all output of preCICE.
   *
   * If done after configuration, this overwrites the output location specified
   * in the configuration.
   */
//  void setExportLocation (
//    const std::string& location,
//    int                exportType = constants::exportAll() );

  /**
   * @brief Writes the contained geometries and spacetree to vtk file.
   *
   * The plotting path has to be specified in the configuration of the
   * accessing participant.
   *
   * @param filenamePrefix [IN] Prefix of all plotted files
   */
  void exportMesh (
    const std::string& filenameSuffix,
    int                exportType = constants::exportAll() );


  /**
   * @brief Scales data values according to configuration.
   *
   * Currently, the only scaling supported is a division of the data values
   * through the surface area belonging to its "support". This allows to come
   * from forces to stresses, e.g..
   */
//  void scaleReadData ()

  /**
   * @brief Returns a handle to a created mesh.
   */
  MeshHandle getMeshHandle ( const std::string& meshName );

  /**
   * @brief Runs the solver interface in server mode.
   */
  void runServer();

private:

  struct M2NWrap {
    m2n::M2N::SharedPointer m2n;
    bool isRequesting;
  };

  // @brief Used for writing debug information.
  static tarch::logging::Log _log;

  std::string _accessorName;

  int _accessorProcessRank;

  int _accessorCommunicatorSize;

  impl::PtrParticipant _accessor;

  // @brief Spatial dimensions of problem.
  int _dimensions;

  // @brief If true, the SolverInterfaceImpl is used by one solver only.
  bool _geometryMode;

  // @brief If true, the simulation continues from a written checkpoint.
  bool _restartMode;

  // @brief If true, the interface is run as server for another interface
  bool _serverMode;

  // @brief If true, the interface uses a server to operate on coupling data.
  bool _clientMode;

  // @brief Communication when for client-server mode.
  //com::Communication::SharedPointer _clientServerCommunication;

  // @brief Geometry name to mesh ID mapping.
  std::map<std::string,int> _meshIDs;

  //@brief dataIDs referenced by meshID and data name
  std::map<int,std::map<std::string,int> > _dataIDs;

  // @brief For plotting of used mesh neighbor-relations
  query::ExportVTKNeighbors _exportVTKNeighbors;

  std::map<std::string,M2NWrap> _m2ns;

  // @brief Holds information about solvers participating in the coupled simulation.
  std::vector<impl::PtrParticipant> _participants;

  cplscheme::PtrCouplingScheme _couplingScheme;

  // @brief To avoid multiple plots during one subcycled global timestep.
  //bool _alreadyPlotThisTimestep;

  // @brief Defines a timestep interval for global checkpointing
  int _checkpointTimestepInterval;

  // @brief Name of checkpoint file.
  std::string _checkpointFileName;

  // @brief Counts calls to advance for plotting.
  long int _numberAdvanceCalls;

//  // @brief Locks the next receive operation of the server to a specific client.
//  int _lockServerToClient;

  // @brief Manages client-server requests, when a server is used.
  RequestManager* _requestManager;

  // @brief In case of a server lock (_lockServerToClient), a specific request
  //        is expected.
  //int _expectRequest;

  void configureM2Ns ( const m2n::M2NConfiguration::SharedPointer& config );

  /**
   * @brief Exports geometries with data and watch point data.
   */
  void handleExports();

  /**
   * @brief Adds exchanged data ids related to accessor to the coupling scheme.
   *
   * From the configuration it is only known which data a participant writes,
   * which data he receives, and which data is exchanged within the coupling
   * scheme. The data actually being sent and received in the coupling scheme,
   * is defined by the intersection of data being written or read by the
   * accessor and the data exchanged in the coupling scheme, respectively.
   *
   * Prerequesits:
   * - _accessor is valid
   * - _couplingScheme is valid
   */
  void addDataToCouplingScheme (
    const cplscheme::CouplingSchemeConfiguration& config );

  /**
   * @brief Configures _couplingScheme to be ready to use.
   *
   * Prerequesits:
   * - _couplingScheme holds a pointer to a valid coupling scheme object
   * - all meshes of geometries are created
   * - all data is added to _data
   * - _accessor points to the accessor
   */
  void addMeshesToCouplingScheme();

  /**
   * @brief Returns true, if the accessor uses the geometry with given name.
   */
  bool isUsingMesh ( const std::string& meshName );

  /**
   * @brief Determines participants providing meshes to other participants.
   */
  void configureSolverGeometries (
    const m2n::M2NConfiguration::SharedPointer& m2nConfig );

  /**
   * @brief Creates the mesh and context data structure of a geometry.
   */
  void createMeshContext ( impl::MeshContext& meshContext );

  /**
   * @brief Computes, performs, and resets all suitable write mappings.
   */
  void mapWrittenData();

  /**
   * @brief Computes, performs, and resets all suitable read mappings.
   */
  void mapReadData();

  /**
   * @brief Performs all data actions with given timing.
   *
   * @param dt [IN] Last timestep length computed by solver.
   * @param partFullDt [IN] Part of current full timestep computed by solver.
   * @param fullDt [IN] Length of current full timestep.
   */
  void performDataActions (
    const std::set<action::Action::Timing>& timings,
    double                 time,
    double                 dt,
    double                 partFullDt,
    double                 fullDt );

  /**
   * @brief Resets written data, displacements and mesh neighbors to export.
   */
  void resetWrittenData();

  /**
   * @brief Sets all dataContext value indices to zero.
   */
  //void resetDataIndices();

  /**
   * @brief Determines participant accessing this interface from the configuration.
   */
  impl::PtrParticipant determineAccessingParticipant (
    const config::SolverInterfaceConfiguration& config );

  /**
   * @brief Selects meshes to be queried from mesh IDs requested.
   *
   * Not all meshes need to be queried separately, if some share the same
   * spacetree for example. This method selects all necessary meshes necessary
   * to be queried.
   */
  void selectInquiryMeshIDs (
      const std::set<int>& meshIDs,
      std::vector<int>&    selectedMeshIDs ) const;

  int markedSkip() const { return 0; }
  int markedQueryDirectly() const { return 1; }
  int markedQuerySpacetree() const { return 2; }

  /**
   * @brief Initializes communication between data server and client.
   */
  void initializeClientServerCommunication();

  /**
   * @brief Initializes communication between master and slaves.
   */
  void initializeMasterSlaveCommunication();

  /**
   * @brief syncs the timestep between slaves and master (all timesteps should be the same!)
   */
  void syncTimestep(double computedTimestepLength);

  // @brief To allow white box tests.
  friend class tests::SolverInterfaceTest;
  friend class tests::SolverInterfaceTestRemote;
  friend class tests::SolverInterfaceTestGeometry;
};

}} // namespace precice, impl

#endif /* PRECICE_IMPL_SOLVERINTERFACEIMPL_HPP_ */
