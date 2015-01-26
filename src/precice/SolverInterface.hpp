// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_SOLVERINTERFACE_HPP_
#define PRECICE_SOLVERINTERFACE_HPP_

#include "MeshHandle.hpp"
#include "Constants.hpp"
#include "ClosestMesh.hpp"
#include "VoxelPosition.hpp"
#include <string>
#include <vector>
#include <set>

/**
 * Pre-declarations.
 */
namespace precice {
  namespace impl {
    class SolverInterfaceImpl;
  }
  namespace tests {
    class SolverInterfaceTest;
    class SolverInterfaceTestRemote;
    class SolverInterfaceTestGeometry;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {

/**
 * @brief Interface to be used by solvers when using preCICE.
 *
 * A solver, i.e. simulation program, using preCICE has to use SolverInterface such
 *
 * -# Create an object of SolverInterface with SolverInterface()
 * -# Configure the object with SolverInterface::configure()
 * -# Initialize preCICE with SolverInterface::initialize()
 * -# Advance to the next (time)step with SolverInterface::advance()
 * -# Finalize preCICE with SolverInterface::finalize()
 */
class SolverInterface
{
public:

  /**
   * @brief Constructor.
   *
   * @param solverName [IN] Name of the solver using the interface. Has to
   *        match the name given for a participant in the xml configuration file.
   * @param solverProcessIndex [IN] If the solver code runs with several processes,
   *        each process using preCICE has specify its index, which has to start
   *        from 0 and end with solverProcessSize - 1.
   * @param solverProcessSize [IN] The number of solver processes using preCICE.
   */
  SolverInterface (
    const std::string& solverName,
    int                solverProcessIndex,
    int                solverProcessSize );

  /**
   * @brief Destructor.
   */
  ~SolverInterface();

  /**
   * @brief Configures preCICE from the given xml file.
   *
   * Only after the configuration a usable state of a SolverInterface
   * object is achieved. However, most of the functionalities in preCICE can be
   * used only after initialize() has been called. Some actions, e.g. specifying
   * the solvers interface mesh, have to be done before initialize is called.
   *
   * In configure, the following is done:
   * - The XML configuration for preCICE is parsed and all objects containing
   *   data are created, but not necessarily filled with data.
   * - If a server is used, communication to that server is established.
   *
   * @param configurationFileName [IN] Name (with path) of the xml configuration
   *        file to be read.
   */
  void configure ( const std::string& configurationFileName );

  /**
   * @brief Fully initializes preCICE to be used.
   *
   * Preconditions:
   * - configure() has been called successfully.
   *
   * Postconditions:
   * - Communication to the coupling partner/s is setup.
   * - Meshes are created from geometries or sent/received to/from coupling partners.
   * - If the solver is not starting the simulation, coupling data is received
   *   from the coupling partner's first computation.
   * - The length limitation of the first solver timestep is computed and returned.
   *
   * @return Maximum length of first timestep to be computed by the solver.
   */
  double initialize();

  /**
   * @brief Initializes coupling data of an implicit coupling scheme.
   *
   * When in a coupled simulation an implicit coupling scheme is used, the
   * starting values for the coupling data are assumed to be zero by default. If
   * this is not the desired behavior, this method can be used to specify values
   * different from zero using the write data methods. Only the first participant
   * of the coupled simulation has to call this method, the second participant
   * receives the values on calling initialize().
   *
   * Preconditions:
   * - initialize() has been called successfully.
   * - The coupling data to be set is written.
   * - advance() has not yet been called.
   * - finalize() has not yet been called.
   *
   * Postconditions:
   * - Initial coupling data values are sent to the coupling partner.
   * - Written coupling data values are reset to zero, in order to allow writing
   *   of values of first computed timestep.
   */
  void initializeData();

  /**
   * @brief Advances preCICE after the solver has computed one timestep.
   *
   * Preconditions:
   * - initialize() has been called successfully.
   * - The solver calling advance() has computed one timestep.
   * - The solver has read and written all coupling data.
   * - The solver has the length of the timestep used available to be passed to
   *   preCICE.
   * - finalize() has not yet been called.
   *
   * Postconditions:
   * - Coupling data values specified to be exchanged in the configuration are
   *   exchanged with coupling partner/s.
   * - Coupling scheme state (computed time, computed timesteps, ...) is updated.
   * - For staggered coupling schemes, the coupling partner has computed one
   *   timestep/iteration with the coupling data written by this participant
   *   available.
   * - The coupling state is printed.
   * - Meshes with data are exported to files if specified in the configuration.
   * - The length limitation of the next solver timestep is computed and returned.
   *
   * @param computedTimestepLength [IN] Length of timestep computed by solver.
   * @return Maximum length of next timestep to be computed by solver.
   */
  double advance ( double computedTimestepLength );

  /**
   * @brief Finalizes preCICE.
   *
   * Preconditions:
   * - initialize() has been called successfully.
   * - isCouplingOngoing() has returned false.
   *
   * Postconditions:
   * - Communication channels are closed.
   * - Meshes and data are deallocated
   */
  void finalize();

  /**
   * @brief Returns the number of spatial dimensions configured.
   *
   * Currently, two and three dimensional problems can be solved using preCICE.
   * The dimension is specified in the XML configuration.
   *
   * Preconditions:
   * - configure() has been called successfully.
   */
  int getDimensions() const;

  /**
   * @brief Returns true, if the coupled simulation is still ongoing.
   *
   * Preconditions:
   * - initialize() has been called successfully.
   */
  bool isCouplingOngoing();

  /**
   * @brief Returns true, if new data to be read is available.
   *
   * Data is classified to be new, if it has been received while calling
   * initialize() and before calling advance(), or in the last call of advance().
   * This is always true, if a participant does not make use of subcycling, i.e.
   * choosing smaller timesteps than the limits returned in intitialize() and
   * advance().
   *
   * Preconditions:
   * - initialize() has been called successfully.
   */
  bool isReadDataAvailable();

  /**
   * @brief Returns true, if new data has to be written before calling advance().
   *
   * This is always true, if a participant does not make use of subcycling, i.e.
   * choosing smaller timesteps than the limits returned in intitialize() and
   * advance().
   *
   * Preconditions:
   * - initialize() has been called successfully.
   */
  bool isWriteDataRequired ( double computedTimestepLength );

  /**
   * @brief Returns true, if the current coupling timestep is completed.
   *
   * The following reasons require several solver time steps per coupling time
   * step:
   * - A solver chooses to perform subcycling.
   * - An implicit coupling timestep iteration is not yet converged.
   *
   * Preconditions:
   * - initialize() has been called successfully.
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
   * @brief Returns true, if the data with given name is used.
   */
  bool hasData ( const std::string& dataName, int meshID ) const;

  /**
   * @brief Returns data id corresponding to the given name (from configuration)
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
    const double*        inquiredPoint,
    const std::set<int>& meshIDs );

  /**
   * @brief Inquires, whether voxel is inside, outside or on geometry surface
   *
   * Objects touching the voxel surface are not considered to be in the voxel.
   *
   * @param inquiryCenter   [IN] Center point of voxel.
   * @param halfLengthVoxel [IN] Half lengths of voxel sides.
   * @param includeBoundaries [IN] If true, touching objects are included.
   * @param storeContent [IN] If true, included objects are stored in VoxelPosition.
   *
   * @return Result of inquiry.
   */
  VoxelPosition inquireVoxelPosition (
    const double*        inquiryCenter,
    const double*        halfLengthVoxel,
    bool                 includeBoundaries,
    //bool                 storeContent,
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
   * @brief Sets position of surface mesh vertex, returns ID.
   */
  int setMeshVertex (
    int           meshID,
    const double* position );

  /**
   * @brief Returns the number of vertices of a mesh.
   */
  int getMeshVertexSize(int meshID);

  /**
   * @brief Sets several spatial vertex positions.
   *
   * Pre-conditions:
   * - A not incremental write-mapping is configured for the mesh with given
   *   meshID.
   * Post-conditions:
   * - If no write mapping is configured, the ids are not changed.
   * - If a (not incremental) write mapping is configured, the ids are filled.
   *
   * @param ids [OUT] IDs for data to be written from given positions.
   */
  void setMeshVertices (
    int     meshID,
    int     size,
    double* positions,
    int*    ids );

  /**
   * @brief Gets spatial vertex positions for given IDs.
   *
   * @param ids [IN] IDs obtained when setting write positions.
   * @param positions [OUT] Positions corresponding to IDs.
   */
  void getMeshVertices (
    int     meshID,
    int     size,
    int*    ids,
    double* positions );

  /**
   * @brief Gets mesh vertex IDs from positions.
   *
   * @param size [IN] Number of positions, ids.
   * @param positions [IN] Positions (x,y,z,x,y,z,...) to find ids for.
   * @param ids [OUT] IDs corresponding to positions.
   */
  void getMeshVertexIDsFromPositions (
    int     meshID,
    int     size,
    double* positions,
    int*    ids );

  /**
   * @brief Sets surface mesh edge from vertex IDs, returns edge ID.
   */
  int setMeshEdge (
    int meshID,
    int firstVertexID,
    int secondVertexID );

  /**
   * @brief Sets surface mesh triangle from edge IDs.
   */
  void setMeshTriangle (
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID );

  /**
   * @brief Sets surface mesh triangle from vertex IDs.
   *
   * This routine is supposed to be used, when no edge information is available
   * per se. Edges are created on the fly within preCICE. This routine is
   * significantly slower than the one using edge IDs, since it needs to check,
   * whether an edge is created already or not.
   */
  void setMeshTriangleWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID );

  /**
   * @brief Sets surface mesh quadrangle from edge IDs.
   */
  void setMeshQuad (
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID,
    int fourthEdgeID );

  /**
   * @brief Sets surface mesh quadrangle from vertex IDs.
   *
   * This routine is supposed to be used, when no edge information is available
   * per se. Edges are created on the fly within preCICE. This routine is
   * significantly slower than the one using edge IDs, since it needs to check,
   * whether an edge is created already or not.
   */
  void setMeshQuadWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID );

  /**
   * @brief Computes and maps all read data mapped to the mesh with given ID.
   *
   */
  void mapReadDataTo ( int toMeshID );

  /**
   * @brief Computes and maps all write data mapped from the mesh with given ID.
   */
  void mapWriteDataFrom ( int fromMeshID );

  /**
   * @brief Writes vector data values given as block.
   *
   * The block must contain the vector values in the following form:
   * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
   * the number of vector values. In 2D, the z-components are removed.
   *
   * @param dataID [IN] ID of the data to be written.
   * @param size [IN] Number of indices, and number of values * dimensions.
   * @param values [IN] Values of the data to be written.
   */
  void writeBlockVectorData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Write vectorial data to the geometry interface
   *
   * The exact mapping and communication must be specified in XYZ.
   *
   * @param dataID       [IN] ID of the data to be written, e.g. 1 = forces
   * @param dataPosition [IN] Position (coordinate, e.g.) of data to be written
   * @param dataValue    [IN] Value of the data to be written
   */
  void writeVectorData (
    int           dataID,
    int           valueIndex,
    const double* value );


  /**
   * @brief Writes scalar data values given as block.
   *
   * @param dataID [IN] ID of the data to be written.
   * @param size [IN] Number of valueIndices, and number of values.
   * @param values [IN] Values of the data to be written.
   */
  void writeBlockScalarData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Write scalar data to the geometry interface
   *
   * The exact mapping and communication must be specified in XYZ.
   *
   * @param dataID       [IN] ID of the data to be written (2 = temperature, e.g.)
   * @param dataPosition [IN] Position (coordinate, e.g.) of data to be written
   * @param dataValue    [IN] Value of the data to be written
   */
  void writeScalarData (
    int    dataID,
    int    valueIndex,
    double value );

  /**
   * @brief Reads vector data values given as block.
   *
   * The block contains the vector values in the following form:
   * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
   * the number of vector values. In 2D, the z-components are removed.
   *
   * @param dataID [IN] ID of the data to be read.
   * @param size [IN] Number of indices, and number of values * dimensions.
   * @param valueIndices [IN] Indices (from setReadPosition()) of data values.
   * @param values [IN] Values of the data to be read.
   */
  void readBlockVectorData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Reads vector data from the coupling mesh.
   *
   * @param dataID       [IN]  ID of the data to be read, e.g. 1 = forces
   * @param dataPosition [IN]  Position (coordinate, e.g.) of data to be read
   * @param dataValue    [OUT] Read data value
   */
  void readVectorData (
    int     dataID,
    int     valueIndex,
    double* value );

  /**
   * @brief Reads scalar data values given as block.
   *
   * @param dataID [IN] ID of the data to be written.
   * @param size [IN] Number of valueIndices, and number of values.
   * @param values [IN] Values of the data to be written.
   */
  void readBlockScalarData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Read scalar data from the geometry interface.
   *
   * The exact mapping and communication must be specified in XYZ.
   *
   * @param dataID       [IN]  ID of the data to be read, e.g. 2 = temperatures
   * @param dataPosition [IN]  Position (coordinate, e.g.) of data to be read
   * @param dataValue    [OUT] Read data value
   */
  void readScalarData (
    int     dataID,
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

private:

  // @brief Pointer to implementation of SolverInterface behavior.
  impl::SolverInterfaceImpl* _impl;

  /**
   * @brief Disable copy construction by making copy constructor private.
   */
  SolverInterface ( const SolverInterface& copy );

  /**
   * @brief Disable assignment construction by making assign. constructor private.
   */
  SolverInterface& operator= ( const SolverInterface& assign );

  // @brief To allow white box tests.
  friend class tests::SolverInterfaceTest;
  friend class tests::SolverInterfaceTestRemote;
  friend class tests::SolverInterfaceTestGeometry;
};

} // namespace precice

#endif /* PRECICE_SOLVERINTERFACE_HPP_ */
