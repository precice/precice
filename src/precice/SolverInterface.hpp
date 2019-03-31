#pragma once

#include "MeshHandle.hpp"
#include "Constants.hpp"
#include <string>
#include <vector>
#include <set>
#include <memory>

/**
 * Pre-declarations.
 */
namespace precice {
  namespace impl {
    class SolverInterfaceImpl;
  }
  namespace testing {
      struct WhiteboxAccessor;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {

/**
 * @brief Main Application Programming Interface of preCICE
 *
 * To adapt a solver to preCICE, follow the following main structure:
 *
 * -# Create an object of SolverInterface with SolverInterface()
 * -# Configure the object with SolverInterface::configure()
 * -# Initialize preCICE with SolverInterface::initialize()
 * -# Advance to the next (time)step with SolverInterface::advance()
 * -# Finalize preCICE with SolverInterface::finalize()
 *
 *  @note
 *  We use solver, simulation code, and participant as synonyms.
 *  The preferred name in the documentation is participant.
 */
class SolverInterface
{
public:

  ///@name Construction and Configuration
  ///@{

  /**
   * @param[in] participantName Name of the participant using the interface. Has to
   *        match the name given for a participant in the xml configuration file.
   * @param[in] solverProcessIndex If the solver code runs with several processes,
   *        each process using preCICE has to specify its index, which has to start
   *        from 0 and end with solverProcessSize - 1.
   * @param[in] solverProcessSize The number of solver processes using preCICE.
   */
  SolverInterface (
    const std::string& participantName,
    int                solverProcessIndex,
    int                solverProcessSize );

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
   * - Communication between master and slaves is established.
   *
   * @pre configure() has not yet been called
   *
   * @param[in] configurationFileName Name (with path) of the xml configuration file to be read.
   */
  void configure ( const std::string& configurationFileName );

  ///@}

  /// @name Steering Methods
  ///@{

  /**
   * @brief Fully initializes preCICE
   *
   * @pre configure() has been called successfully.
   * @pre initialize() has not yet bee called.
   *
   * @post Parallel communication to the coupling partner/s is setup.
   * @post Meshes are exchanged between coupling partners and the parallel partitions are created.
   * @post [Serial Coupling Scheme] If the solver is not starting the simulation, coupling data is received
   * from the coupling partner's first computation.
   *
   * @return Maximum length of first timestep to be computed by the solver.
   */
  double initialize();

  /**
   * @brief Initializes coupling data.
   *
   * The starting values for coupling data are zero by default.
   *
   * To provide custom values, first set the data using the Data Access methods and
   * call this method to finally exchange the data.
   *
   * \par Serial Coupling Scheme
   * Only the first participant has to call this method, the second participant
   * receives the values on calling initialize().
   *
   * \par Parallel Coupling Scheme
   * Values in both directions are exchanged.
   * Both participants need to call initializeData().
   *
   * @pre initialize() has been called successfully.
   * @pre The action WriteInitialData is required 
   * @pre advance() has not yet been called.
   * @pre finalize() has not yet been called.
   *
   * @post Initial coupling data was exchanged.
   *
   * @see isActionRequired  
   * @see precice::constants::actionWriteInitialData
   */
  void initializeData();

  /**
   * @brief Advances preCICE after the solver has computed one timestep.
   *
   * @param[in] computedTimestepLength Length of timestep used by the solver.
   *
   * @pre initialize() has been called successfully.
   * @pre The solver has computed one timestep.
   * @pre The solver has written all coupling data.
   * @pre finalize() has not yet been called.
   *
   * @post Coupling data values specified in the configuration are exchanged.
   * @post Coupling scheme state (computed time, computed timesteps, ...) is updated.
   * @post The coupling state is logged.
   * @post Configured data mapping schemes are applied.
   * @post [Second Participant] Configured post processing schemes are applied.
   * @post Meshes with data are exported to files if configured.
   *
   * @return Maximum length of next timestep to be computed by solver.
   */
  double advance ( double computedTimestepLength );

  /**
   * @brief Finalizes preCICE.
   *
   * @pre initialize() has been called successfully.
   *
   * @post Communication channels are closed.
   * @post Meshes and data are deallocated
   *
   * @see isCouplingOngoing()
   */
  void finalize();

  ///@}
  
  ///@name Status Queries
  ///@{

  /**
   * @brief Returns the number of spatial dimensions configured.
   *
   * @returns the configured dimension
   *
   * Currently, two and three dimensional problems can be solved using preCICE.
   * The dimension is specified in the XML configuration.
   *
   * @pre configure() has been called successfully.
   */
  int getDimensions() const;

  /**
   * @brief Checks if the coupled simulation is still ongoing.
   *
   * @returns whether the coupling is ongoing.
   *
   * A coupling is ongoing as long as
   * - the maximum number of timesteps has not been reached, and
   * - the final time has not been reached.
   *
   * @pre initialize() has been called successfully.
   *
   * @see advance()
   *
   * @note
   * The user should call finalize() after this function returns false.
   */
  bool isCouplingOngoing();

  /**
   * @brief Checks if new data to be read is available.
   *
   * @returns whether new data is available to be read.
   *
   * Data is classified to be new, if it has been received while calling
   * initialize() and before calling advance(), or in the last call of advance().
   * This is always true, if a participant does not make use of subcycling, i.e.
   * choosing smaller timesteps than the limits returned in intitialize() and
   * advance().
   *
   * @pre initialize() has been called successfully.
   *
   * @note
   * It is allowed to read data even if this function returns false.
   * This is not recommended due to performance reasons.
   * Use this function to prevent unnecessary reads.
   */
  bool isReadDataAvailable();

  /**
   * @brief Checks if new data has to be written before calling advance().
   *
   * @param[in] computedTimestepLength Length of timestep used by the solver.
   *
   * @return whether new data has to be written.
   *
   * This is always true, if a participant does not make use of subcycling, i.e.
   * choosing smaller timesteps than the limits returned in intitialize() and
   * advance().
   *
   * @pre initialize() has been called successfully.
   *
   * @note
   * It is allowed to write data even if this function returns false.
   * This is not recommended due to performance reasons.
   * Use this function to prevent unnecessary writes.
   */
  bool isWriteDataRequired ( double computedTimestepLength );

  /**
   * @brief Checks if the current coupling timestep is completed.
   *
   * @returns whether the timestep is complete.
   *
   * The following reasons require several solver time steps per coupling time
   * step:
   * - A solver chooses to perform subcycling.
   * - An implicit coupling timestep iteration is not yet converged.
   *
   * @pre initialize() has been called successfully.
   */
  bool isTimestepComplete();

  /**
   * @brief Returns whether the solver has to evaluate the surrogate model representation.
   *
   * @deprecated
   * Only necessary for deprecated manifold mapping.
   *
   * @returns whether the surrogate model has to be evaluated.
   *
   * @note
   * The solver may still have to evaluate the fine model representation.
   *
   * @see hasToEvaluateFineModel()
   */
  bool hasToEvaluateSurrogateModel();

  /**
   * @brief Checks if the solver has to evaluate the fine model representation.
   *
   * @deprecated
   * Only necessary for deprecated manifold mapping.
   *
   * @returns whether the fine model has to be evaluated.
   *
   * @note
   * The solver may still have to evaluate the surrogate model representation.
   *
   * @see hasToEvaluateSurrogateModel()
   */
  bool hasToEvaluateFineModel();

  ///@}

  ///@name Action Methods
  ///@{

  /**
   * @brief Checks if the provided action is required.
   *
   * @param[in] action the name of the action
   * @returns whether the action is required
   *
   * Some features of preCICE require a solver to perform specific actions, in
   * order to be in valid state for a coupled simulation. A solver is made
   * eligible to use those features, by querying for the required actions,
   * performing them on demand, and calling fulfilledAction() to signalize
   * preCICE the correct behavior of the solver.
   *
   * @see fulfilledAction()
   * @see cplscheme::constants
   */
  bool isActionRequired ( const std::string& action );

  /**
   * @brief Indicates preCICE that a required action has been fulfilled by a solver.
   *
   * @pre The solver fulfilled the specified action.
   *
   * @param[in] action the name of the action
   *
   * @see requireAction()
   * @see cplscheme::constants
   */
  void fulfilledAction ( const std::string& action );

  ///@}

  ///@name Mesh Access
  ///@anchor precice-mesh-access
  ///@{

  /*
   * @brief Resets mesh with given ID.
   *
   * Has to be called, everytime the positions for data to be mapped
   * changes. Only has an effect, if the mapping used is non-stationary and
   * non-incremental.
   */
//  void resetMesh ( int meshID );

  /**
   * @brief Checks if the mesh with given name is used by a solver.
   *
   * @param[in] meshName the name of the mesh
   * @returns whether the mesh is used.
   */
  bool hasMesh ( const std::string& meshName ) const;

  /**
   * @brief Returns the ID belonging to the mesh with given name.
   * 
   * @param[in] meshName the name of the mesh
   * @returns the id of the corresponding mesh
   */
  int getMeshID ( const std::string& meshName );

  /**
   * @brief Returns a id-set of all used meshes by this participant.
   *
   * @returns the set of ids.
   */
  std::set<int> getMeshIDs();

  /**
   * @brief Returns a handle to a created mesh.
   * 
   * @param[in] meshName the name of the mesh
   * @returns the handle to the mesh
   *
   * @see precice::MeshHandle
   */
  MeshHandle getMeshHandle ( const std::string& meshName );

  /**
   * @brief Creates a mesh vertex
   *
   * @param[in] meshID the id of the mesh to add the vertex to.
   * @param[in] position a pointer to the coordinates of the vertex.
   * @returns the id of the created vertex
   *
   * @pre initialize() has not yet been called
   * @pre count of available elements at position matches the configured dimension
   *
   * @see getDimensions()
   */
  int setMeshVertex (
    int           meshID,
    const double* position );

  /**
   * @brief Returns the number of vertices of a mesh.
   *
   * @param[in] meshID the id of the mesh
   * @returns the amount of the vertices of the mesh
   */
  int getMeshVertexSize(int meshID);

  /**
   * @brief Creates multiple mesh vertices
   *
   * @param[in] meshID the id of the mesh to add the vertices to.
   * @param[in] size Number of vertices to create
   * @param[in] positions a pointer to the coordinates of the vertices
   *            The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   *            The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @param[out] ids The ids of the created vertices
   *
   * @pre initialize() has not yet been called
   * @pre count of available elements at positions matches the configured dimension * size
   * @pre count of available elements at ids matches size
   *
   * @see getDimensions()
   */
  void setMeshVertices (
    int     meshID,
    int     size,
    double* positions,
    int*    ids );

  /**
   * @brief Get vertex positions for multiple vertex ids from a given mesh
   *
   * @param[in] meshID the id of the mesh to read the vertices from.
   * @param[in] size Number of vertices to lookup
   * @param[in] ids The ids of the vertices to lookup
   * @param[out] positions a pointer to memory to write the coordinates to
   *            The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   *            The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @pre count of available elements at positions matches the configured dimension * size
   * @pre count of available elements at ids matches size
   *
   * @see getDimensions()
   */
  void getMeshVertices (
    int     meshID,
    int     size,
    int*    ids,
    double* positions );

  /**
   * @brief Gets mesh vertex IDs from positions.
   *
   * @param[in] meshID ID of the mesh to retrieve positions from
   * @param[in] size Number of vertices to lookup.
   * @param[in] positions Positions to find ids for.
   *            The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   *            The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   * @param[out] ids IDs corresponding to positions.
   *
   * @pre count of available elements at positions matches the configured dimension * size
   * @pre count of available elements at ids matches size
   *
   * @note prefer to reuse the IDs returned from calls to setMeshVertex() and setMeshVertices().
   */
  void getMeshVertexIDsFromPositions (
    int     meshID,
    int     size,
    double* positions,
    int*    ids );

  /**
   * @brief Sets mesh edge from vertex IDs, returns edge ID.
   *
   * @param[in] meshID ID of the mesh to add the edge to
   * @param[in] firstVertexID ID of the first vertex of the edge
   * @param[in] secondVertexID ID of the second vertex of the edge
   *
   * @return the ID of the edge
   *
   * @pre vertices with firstVertexID and secondVertexID were added to the mesh with the ID meshID
   */
  int setMeshEdge (
    int meshID,
    int firstVertexID,
    int secondVertexID );

  /**
   * @brief Sets mesh triangle from edge IDs
   *
   * @param[in] meshID ID of the mesh to add the triangle to
   * @param[in] firstEdgeID ID of the first edge of the triangle
   * @param[in] secondEdgeID ID of the second edge of the triangle
   * @param[in] thirdEdgeID ID of the third edge of the triangle
   *
   * @pre edges with firstEdgeID, secondEdgeID, and thirdEdgeID were added to the mesh with the ID meshID
   */
  void setMeshTriangle (
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID );

  /**
   * @brief Sets mesh triangle from vertex IDs.
   *
   * @warning
   * This routine is supposed to be used, when no edge information is available
   * per se. Edges are created on the fly within preCICE. This routine is
   * significantly slower than the one using edge IDs, since it needs to check,
   * whether an edge is created already or not.
   * 
   * @param[in] meshID ID of the mesh to add the triangle to
   * @param[in] firstVertexID ID of the first vertex of the triangle
   * @param[in] secondVertexID ID of the second vertex of the triangle
   * @param[in] thirdVertexID ID of the third vertex of the triangle
   *
   * @pre edges with firstVertexID, secondVertexID, and thirdVertexID were added to the mesh with the ID meshID
   */
  void setMeshTriangleWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID );

  /**
   * @brief Sets mesh Quad from edge IDs.
   *
   * @param[in] meshID ID of the mesh to add the Quad to
   * @param[in] firstEdgeID ID of the first edge of the Quad
   * @param[in] secondEdgeID ID of the second edge of the Quad
   * @param[in] thirdEdgeID ID of the third edge of the Quad
   * @param[in] fourthEdgeID ID of the forth edge of the Quad
   *
   * @pre edges with firstEdgeID, secondEdgeID, thirdEdgeID, and fourthEdgeID were added to the mesh with the ID meshID
   *
   * @warning Quads are not fully implemented yet.
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
   * @warning
   * This routine is supposed to be used, when no edge information is available
   * per se. Edges are created on the fly within preCICE. This routine is
   * significantly slower than the one using edge IDs, since it needs to check,
   * whether an edge is created already or not.
   * 
   * @param[in] meshID ID of the mesh to add the Quad to
   * @param[in] firstVertexID ID of the first vertex of the Quad
   * @param[in] secondVertexID ID of the second vertex of the Quad
   * @param[in] thirdVertexID ID of the third vertex of the Quad
   * @param[in] fourthVertexID ID of the fourth vertex of the Quad
   *
   * @pre edges with firstVertexID, secondVertexID, thirdVertexID, and fourthVertexID were added to the mesh with the ID meshID
   *
   */
  void setMeshQuadWithEdges (
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID );

  ///@}

  ///@name Data Access
  ///@{

  /**
   * @brief Checks if the data with given name is used by a solver and mesh.
   *
   * @param[in] dataName the name of the data
   * @param[in] meshID the id of the associated mesh
   * @returns whether the mesh is used.
   */
  bool hasData ( const std::string& dataName, int meshID ) const;

  /**
   * @brief Returns the ID of the data associated with the given name and mesh.
   * 
   * @param[in] dataName the name of the data
   * @param[in] meshID the id of the associated mesh
   *
   * @returns the id of the corresponding data
   */
  int getDataID ( const std::string& dataName, int meshID );

  /**
   * @brief Computes and maps all read data mapped to the mesh with given ID.
   *
   * This is an explicit request to map read data to the Mesh associated with toMeshID.
   * It also computes the mapping if necessary.
   *
   * @pre A mapping to toMeshID was configured.
   */
  void mapReadDataTo ( int toMeshID );

  /**
   * @brief Computes and maps all write data mapped from the mesh with given ID.
   *
   * This is an explicit request to map write data from the Mesh associated with fromMeshID.
   * It also computes the mapping if necessary.
   *
   * @pre A mapping from fromMeshID was configured.
   */
  void mapWriteDataFrom ( int fromMeshID );

  /**
   * @brief Writes vector data given as block.
   *
   * This function writes values of specified vertices to a dataID.
   * Values are provided as a block of continuous memory.
   * valueIndices contains the indices of the vertices
   *
   * The 2D-format of values is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   * The 3D-format of values is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @param[in] dataID ID to write to.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[in] values pointer to the vector values.
   *
   * @pre count of available elements at values matches the configured dimension * size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeBlockVectorData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Writes vector data to a vertex
   *
   * This function writes a value of a specified vertex to a dataID.
   * Values are provided as a block of continuous memory.
   *
   * The 2D-format of value is (x, y)
   * The 3D-format of value is (x, y, z)
   *
   * @param[in] dataID ID to write to.
   * @param[in] valueIndex Index of the vertex.
   * @param[in] value pointer to the vector value.
   *
   * @pre count of available elements at value matches the configured dimension
   * @pre initialize() has been called
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeVectorData (
    int           dataID,
    int           valueIndex,
    const double* value );


  /**
   * @brief Writes scalar data given as block.
   *
   * This function writes values of specified vertices to a dataID.
   * Values are provided as a block of continuous memory.
   * valueIndices contains the indices of the vertices
   *
   * @param[in] dataID ID to write to.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[in] values pointer to the values.
   *
   * @pre count of available elements at values matches the given size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeBlockScalarData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Writes scalar data to a vertex
   *
   * This function writes a value of a specified vertex to a dataID.
   *
   * @param[in] dataID ID to write to.
   * @param[in] valueIndex Index of the vertex.
   * @param[in] value the value to write.
   *
   * @pre initialize() has been called
   *
   * @see SolverInterface::setMeshVertex()
   */
  void writeScalarData (
    int    dataID,
    int    valueIndex,
    double value );

  /**
   * @brief Reads vector data into a provided block.
   *
   * This function reads values of specified vertices from a dataID.
   * Values are read into a block of continuous memory.
   * valueIndices contains the indices of the vertices.
   *
   * The 2D-format of values is (d0x, d0y, d1x, d1y, ..., dnx, dny)
   * The 3D-format of values is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
   *
   * @param[in] dataID ID to read from.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[out] values pointer to read destination.
   *
   * @pre count of available elements at values matches the configured dimension * size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @post values contain the read data as specified in the above format.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readBlockVectorData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Reads vector data form a vertex
   *
   * This function reads a value of a specified vertex from a dataID.
   * Values are provided as a block of continuous memory.
   *
   * The 2D-format of value is (x, y)
   * The 3D-format of value is (x, y, z)
   *
   * @param[in] dataID ID to read from.
   * @param[in] valueIndex Index of the vertex.
   * @param[out] value pointer to the vector value.
   *
   * @pre count of available elements at value matches the configured dimension
   * @pre initialize() has been called
   *
   * @post value contains the read data as specified in the above format.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readVectorData (
    int     dataID,
    int     valueIndex,
    double* value );

  /**
   * @brief Reads scalar data as a block.
   *
   * This function reads values of specified vertices from a dataID.
   * Values are provided as a block of continuous memory.
   * valueIndices contains the indices of the vertices.
   *
   * @param[in] dataID ID to read from.
   * @param[in] size Number n of vertices.
   * @param[in] valueIndices Indices of the vertices.
   * @param[out] values pointer to the read destination.
   *
   * @pre count of available elements at values matches the given size
   * @pre count of available elements at valueIndices matches the given size
   * @pre initialize() has been called
   *
   * @post values contains the read data.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readBlockScalarData (
    int     dataID,
    int     size,
    int*    valueIndices,
    double* values );

  /**
   * @brief Reads scalar data of a vertex.
   *
   * This function reads a value of a specified vertex from a dataID.
   *
   * @param[in] dataID ID to read from.
   * @param[in] valueIndex Index of the vertex.
   * @param[out] value read destination of the value.
   *
   * @pre initialize() has been called
   *
   * @post value contains the read data.
   *
   * @see SolverInterface::setMeshVertex()
   */
  void readScalarData (
    int     dataID,
    int     valueIndex,
    double& value );

  ///@}

private:

  /// Pointer to implementation of SolverInterface.
  std::unique_ptr<impl::SolverInterfaceImpl> _impl;

  /// Disable copy construction by making copy constructor private.
  SolverInterface ( const SolverInterface& copy );

  /// Disable assignment construction by making assign. constructor private.
  SolverInterface& operator= ( const SolverInterface& assign );

  // @brief To allow white box tests.
  friend struct testing::WhiteboxAccessor;
};

} // namespace precice
