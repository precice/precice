#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief C language bindings to main Application Programming Interface of preCICE
 *
 */

///@name Construction and Configuration
///@{

/**
 * @param[in] participantName Name of the participant using the interface. Has to
 *        match the name given for a participant in the xml configuration file.
 * @param[in] configurationFileName Name (with path) of the xml configuration file.
 * @param[in] solverProcessIndex If the solver code runs with several processes,
 *        each process using preCICE has to specify its index, which has to start
 *        from 0 and end with solverProcessSize - 1.
 * @param[in] solverProcessSize The number of solver processes using preCICE.
 * @param[in] communicator A pointer to an MPI_Comm to use as communicator.
 */
void precicec_createSolverInterface_withCommunicator(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize,
    void *      communicator);

/**
 * @brief Creates the coupling interface and configures it.
 *
 * Has to be called before any other method of this interface.
 *
 * @param[in] participantName Name of the participant accessing the interface. Has to
 *                          match one of the names specified in the
 *                          configuration xml file.
 * @param[in] configFileName (Path and) name of the xml configuration file
 *                            containing the precice configuration.
 * @param[in] solverProcessIndex If the solver code runs with several processes,
 *                               each process using preCICE has to specify its index, which has to start
 *                               from 0 and end with solverProcessSize - 1.
 * @param[in] solverProcessSize The number of solver processes using preCICE.
 */
void precicec_createSolverInterface(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize);

///@}

/// @name Steering Methods
///@{

/**
 * @brief Initiates the coupling to the coupling supervisor.
 *
 * @return Maximal length of first timestep to be computed by solver.
 */
double precicec_initialize();

/**
 * @brief Initializes coupling data.
 */
void precicec_initialize_data();

/**
 * @brief Exchanges data between solver and coupling supervisor.
 *
 * @param[in] computedTimestepLength Length of timestep computed by solver.
 * @return Maximal length of next timestep to be computed by solver.
 */
double precicec_advance(double computedTimestepLength);

/**
 * @brief Finalizes the coupling to the coupling supervisor.
 */
void precicec_finalize();

///@}

///@name Status Queries
///@{

/**
 * @brief Returns the number of spatial configurations for the coupling.
 */
int precicec_getDimensions();

/**
 * @brief Returns true (->1), if the coupled simulation is ongoing
 */
int precicec_isCouplingOngoing();

/**
 * @brief Returns true (->1), if new data to read is available.
 */
int precicec_isReadDataAvailable();

/**
 * @brief Checks if new data has to be written before calling advance().
 *
 * @param[in] computedTimestepLength Length of timestep used by the solver.
 *
 * @return true (->1) if new data has to be written.
 */
int precicec_isWriteDataRequired(double computedTimestepLength);

/**
 * @brief Returns true (->1), if the coupling time window is completed.
 */
int precicec_isTimeWindowComplete();

/**
 * @brief Returns whether the solver has to evaluate the surrogate model representation.
 */
int precicec_hasToEvaluateSurrogateModel();

/**
 * @brief Returns whether the solver has to evaluate the fine model representation.
 */
int precicec_hasToEvaluateFineModel();

///@}

///@name Action Methods
///@{

/**
 * @brief Checks if the provided action is required.
 * @param[in] action the name of the action
 * @returns whether the action is required
 */
int precicec_isActionRequired(const char *action);

/**
 * @brief Indicates preCICE that a required action has been fulfilled by a solver.
 * @pre The solver fulfilled the specified action.
 *
 * @param[in] action the name of the action
 */
void precicec_markActionFulfilled(const char *action);

///@}

///@name Mesh Access
///@anchor precice-mesh-access
///@{

/**
 * @brief Checks if the mesh with given name is used by a solver.
 *
 * @param[in] meshName the name of the mesh
 * @returns whether the mesh is used.
 */
int precicec_hasMesh(const char *meshName);

/**
 * @brief Returns id belonging to the given mesh name
 */
int precicec_getMeshID(const char *meshName);

/**
 * @brief Creates a mesh vertex
 *
 * @param[in] meshID the id of the mesh to add the vertex to.
 * @param[in] position a pointer to the coordinates of the vertex.
 * @returns the id of the created vertex
 */
int precicec_setMeshVertex(
    int           meshID,
    const double *position);

/**
 * @brief Returns the number of vertices of a mesh.
 *
 * @param[in] meshID the id of the mesh
 * @returns the amount of the vertices of the mesh
 */
int precicec_getMeshVertexSize(int meshID);

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
 */
void precicec_setMeshVertices(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids);

/**
 * @brief Get vertex positions for multiple vertex ids from a given mesh
 *
 * @param[in] meshID the id of the mesh to read the vertices from.
 * @param[in] size Number of vertices to lookup
 * @param[in] ids The ids of the vertices to lookup
 * @param[out] positions a pointer to memory to write the coordinates to
 *            The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
 *            The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
 */
void precicec_getMeshVertices(
    int        meshID,
    int        size,
    const int *ids,
    double *   positions);

/**
 * @brief Gets mesh vertex IDs from positions.
 *
 * @param[in] meshID ID of the mesh to retrieve positions from
 * @param[in] size Number of vertices to lookup.
 * @param[in] positions Positions to find ids for.
 *            The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
 *            The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
 * @param[out] ids IDs corresponding to positions.
 */
void precicec_getMeshVertexIDsFromPositions(
    int           meshID,
    int           size,
    const double *positions,
    int *         ids);

/**
 * @brief Sets mesh edge from vertex IDs, returns edge ID.
 *
 * @param[in] meshID ID of the mesh to add the edge to
 * @param[in] firstVertexID ID of the first vertex of the edge
 * @param[in] secondVertexID ID of the second vertex of the edge
 *
 * @return the ID of the edge
 */
int precicec_setMeshEdge(
    int meshID,
    int firstVertexID,
    int secondVertexID);

/**
 * @brief Sets mesh triangle from edge IDs
 *
 * @param[in] meshID ID of the mesh to add the triangle to
 * @param[in] firstEdgeID ID of the first edge of the triangle
 * @param[in] secondEdgeID ID of the second edge of the triangle
 * @param[in] thirdEdgeID ID of the third edge of the triangle
 */
void precicec_setMeshTriangle(
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID);

/**
 * @brief Sets a triangle from vertex IDs. Creates missing edges.
 */
void precicec_setMeshTriangleWithEdges(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID);

/**
 * @brief Sets mesh Quad from edge IDs.
 *
 * @param[in] meshID ID of the mesh to add the Quad to
 * @param[in] firstEdgeID ID of the first edge of the Quad
 * @param[in] secondEdgeID ID of the second edge of the Quad
 * @param[in] thirdEdgeID ID of the third edge of the Quad
 * @param[in] fourthEdgeID ID of the forth edge of the Quad
 */
void precicec_setMeshQuad(
    int meshID,
    int firstEdgeID,
    int secondEdgeID,
    int thirdEdgeID,
    int fourthEdgeID);

/**
  * @brief Sets surface mesh quadrangle from vertex IDs.
  *
  * @param[in] meshID ID of the mesh to add the Quad to
  * @param[in] firstVertexID ID of the first vertex of the Quad
  * @param[in] secondVertexID ID of the second vertex of the Quad
  * @param[in] thirdVertexID ID of the third vertex of the Quad
  * @param[in] fourthVertexID ID of the fourth vertex of the Quad
 */
void precicec_setMeshQuadWithEdges(
    int meshID,
    int firstVertexID,
    int secondVertexID,
    int thirdVertexID,
    int fourthVertexID);

///@}

///@name Data Access
///@{

/**
 * @brief Returns true (!=0), if data with given name is available.
 */
int precicec_hasData(const char *dataName, int meshID);

/**
 * @brief Returns the data id belonging to the given name.
 *
 * The given name (dataName) has to be one of the names specified in the
 * configuration file. The data id obtained can be used to read and write
 * data to and from the coupling mesh.
 */
int precicec_getDataID(const char *dataName, int meshID);

/**
 * @brief Computes and maps all read data mapped to mesh with given ID.
 */
void precicec_mapReadDataTo(int toMeshID);

/**
 * @brief Computes and maps all write data mapped from mesh with given ID.
 */
void precicec_mapWriteDataFrom(int fromMeshID);

/**
 * @brief Writes vector data values given as block.
 *
 * The block must contain the vector values in the following form:
 * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
 * the number of vector values. In 2D, the z-components are removed.
 *
 * @param[in] dataID ID of the data to be written.
 * @param[in] size Number of indices, and number of values * dimensions.
 * @param[in] values Values of the data to be written.
 */
void precicec_writeBlockVectorData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values);

/**
 * @brief Writes vectorial foating point data to the coupling mesh.
 *
 * @param[in] dataID ID of the data to be written. Obtained by getDataID().
 * @param[in] dataPosition Spatial position of the data to be written.
 * @param[in] dataValue Vectorial data value to be written.
 */
void precicec_writeVectorData(
    int           dataID,
    int           valueIndex,
    const double *dataValue);

/**
 * @brief See precice::SolverInterface::writeBlockScalarData().
 */
void precicec_writeBlockScalarData(
    int           dataID,
    int           size,
    const int *   valueIndices,
    const double *values);

/**
 * @brief Writes scalar floating point data to the coupling mesh.
 *
 * @param[in] dataID ID of the data to be written. Obtained by getDataID().
 * @param[in] dataPosition Spatial position of the data to be written.
 * @param[in] dataValue Scalar data value to be written.
 */
void precicec_writeScalarData(
    int    dataID,
    int    valueIndex,
    double dataValue);

/**
 * @brief Reads vector data values given as block.
 *
 * The block contains the vector values in the following form:
 * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
 * the number of vector values. In 2D, the z-components are removed.
 *
 * @param[in] dataID ID of the data to be read.
 * @param[in] size  Number of indices, and number of values * dimensions.
 * @param[in] valueIndices Indices (from setReadPosition()) of data values.
 * @param[in] values Values of the data to be read.
 */
void precicec_readBlockVectorData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values);

/**
 * @brief Reads vectorial foating point data from the coupling mesh.
 *
 * @param[in] dataID ID of the data to be read. Obtained by getDataID().
 * @param[in] dataPosition Position where the read data should be mapped to.
 * @param[out] dataValue Vectorial data value read.
 */
void precicec_readVectorData(
    int     dataID,
    int     valueIndex,
    double *dataValue);

/**
 * @brief See precice::SolverInterface::readBlockScalarData().
 */
void precicec_readBlockScalarData(
    int        dataID,
    int        size,
    const int *valueIndices,
    double *   values);

/**
 * @brief Reads scalar foating point data from the coupling mesh.
 *
 * @param[in] dataID ID of the data to be read. Obtained by getDataID().
 * @param[in] dataPosition Position where the read data should be mapped to.
 * @param[out] dataValue Scalar data value read.
 */
void precicec_readScalarData(
    int     dataID,
    int     valueIndex,
    double *dataValue);

/** 
 * @brief Returns information on the version of preCICE.
 *
 * Returns a semicolon-separated C-string containing:
 * 
 * 1) the version of preCICE
 * 2) the revision information of preCICE
 * 3) the configuration of preCICE including MPI, PETSC, PYTHON
 */
const char *precicec_getVersionInformation();

// @brief Name of action for writing initial data.
const char *precicec_actionWriteInitialData();

// @brief Name of action for writing iteration checkpoint
const char *precicec_actionWriteIterationCheckpoint();

// @brief Name of action for reading iteration checkpoint.
const char *precicec_actionReadIterationCheckpoint();

#ifdef __cplusplus
}
#endif
