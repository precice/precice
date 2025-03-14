#pragma once

#include <precice/Version.h>
#include <precice/export.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file
 * @brief C language bindings to main Application Programming Interface of preCICE
 *
 * The C bindings are thin wrappers around the C++ API.
 *
 * Refer to \ref precice::Participant for detailed documentation.
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
PRECICE_API void precicec_createParticipant_withCommunicator(
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
PRECICE_API void precicec_createParticipant(
    const char *participantName,
    const char *configFileName,
    int         solverProcessIndex,
    int         solverProcessSize);

///@}

/// @name Steering Methods
///@{

/**
 * @brief Initiates the coupling to the coupling supervisor and initializes coupling data.
 */
PRECICE_API void precicec_initialize();

/**
 * @brief Exchanges data between solver and coupling supervisor.
 *
 * @param[in] computedTimeStepSize Size of time step computed by solver.
 */
PRECICE_API void precicec_advance(double computedTimeStepSize);

/**
 * @brief Finalizes the coupling to the coupling supervisor.
 */
PRECICE_API void precicec_finalize();

///@}

///@name Implicit coupling
///@{

/// @copydoc precice::Participant::requiresWritingCheckpoint()
PRECICE_API int precicec_requiresWritingCheckpoint();

/// @copydoc precice::Participant::requiresReadingCheckpoint()
PRECICE_API int precicec_requiresReadingCheckpoint();

///@}

///@name Status Queries
///@{

/// @copydoc precice::Participant::getMeshDimensions()
PRECICE_API int precicec_getMeshDimensions(const char *meshName);

/// @copydoc precice::Participant::getDataDimensions()
PRECICE_API int precicec_getDataDimensions(const char *meshName, const char *dataName);

/**
 * @brief Returns true (->1), if the coupled simulation is ongoing
 */
PRECICE_API int precicec_isCouplingOngoing();

/**
 * @brief Returns true (->1), if the coupling time window is completed.
 */
PRECICE_API int precicec_isTimeWindowComplete();

/**
 * @brief Returns maximum allowed time step size
 */
PRECICE_API double precicec_getMaxTimeStepSize();

///@}

///@name Mesh Access
///@anchor precice-mesh-access
///@{

/// @copydoc precice::Participant::requiresMeshConnectivityFor()
PRECICE_API int precicec_requiresMeshConnectivityFor(const char *meshName);

/**
 * @brief Creates a mesh vertex
 *
 * @param[in] meshName the name of the mesh to add the vertex to.
 * @param[in] coordinates a pointer to the coordinates of the vertex.
 * @returns the id of the created vertex
 */
PRECICE_API int precicec_setMeshVertex(
    const char *  meshName,
    const double *coordinates);

/**
 * @brief Returns the number of vertices of a mesh.
 *
 * @param[in] meshName the name of the mesh.
 * @returns the amount of the vertices of the mesh
 */
PRECICE_API int precicec_getMeshVertexSize(const char *meshName);

/**
 * @brief Creates multiple mesh vertices
 *
 * @param[in] meshName the name of the mesh to add the vertices to.
 * @param[in] size Number of vertices to create
 * @param[in] coordinates a pointer to the coordinates of the vertices
 *            The 2D-format is (d0x, d0y, d1x, d1y, ..., dnx, dny)
 *            The 3D-format is (d0x, d0y, d0z, d1x, d1y, d1z, ..., dnx, dny, dnz)
 *
 * @param[out] ids The ids of the created vertices
 */
PRECICE_API void precicec_setMeshVertices(
    const char *  meshName,
    int           size,
    const double *coordinates,
    int *         ids);

/**
 * @brief Sets mesh edge from vertex IDs
 *
 * @param[in] meshName the name of the mesh to add the edge to
 * @param[in] firstVertexID ID of the first vertex of the edge
 * @param[in] secondVertexID ID of the second vertex of the edge
 *
 * @pre ids were added to the mesh with the name \p meshName
 */
PRECICE_API void precicec_setMeshEdge(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID);

/**
 * @brief Sets multiple mesh edge from vertex IDs
 *
 * @param[in] meshName the name of the mesh to add the edges to
 * @param[in] size the amount of edges to set
 * @param[in] ids an array containing 2*size vertex IDs
 *
 * @pre ids were added to the mesh with the name \p meshName
 */
PRECICE_API void precicec_setMeshEdges(
    const char *meshName,
    int         size,
    const int * ids);

/**
 * @brief Sets a triangle from vertex IDs. Creates missing edges.
 */
PRECICE_API void precicec_setMeshTriangle(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID);

/**
 * @brief Sets multiple mesh triangles from vertex IDs
 *
 * @param[in] meshName the name of the mesh to add the triangles to
 * @param[in] size the amount of triangles to set
 * @param[in] ids an array containing 3*size vertex IDs
 *
 * @pre vertices in \p ids were added to the mesh with the name \p meshName
 */
PRECICE_API void precicec_setMeshTriangles(
    const char *meshName,
    int         size,
    const int * ids);

/**
 * @brief Sets surface mesh quadrangle from vertex IDs.
 *
 * @param[in] meshName the name of the mesh to add the Quad to
 * @param[in] firstVertexID ID of the first vertex of the Quad
 * @param[in] secondVertexID ID of the second vertex of the Quad
 * @param[in] thirdVertexID ID of the third vertex of the Quad
 * @param[in] fourthVertexID ID of the fourth vertex of the Quad
 */
PRECICE_API void precicec_setMeshQuad(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID,
    int         fourthVertexID);

/**
 * @brief Sets multiple mesh quads from vertex IDs
 *
 * @param[in] meshName the name of the mesh to add the quad to
 * @param[in] size the amount of quads to set
 * @param[in] ids an array containing 4*size vertex IDs
 *
 * @pre vertices in \ids were added to the mesh with the name \p meshName
 */
PRECICE_API void precicec_setMeshQuads(
    const char *meshName,
    int         size,
    const int * ids);

/**
 * @brief Sets mesh tetrahedron from vertex IDs.
 *
 * @param[in] meshName the name of the mesh to add the Tetra to
 * @param[in] firstVertexID ID of the first vertex of the Tetra
 * @param[in] secondVertexID ID of the second vertex of the Tetra
 * @param[in] thirdVertexID ID of the third vertex of the Tetra
 * @param[in] fourthVertexID ID of the fourth vertex of the Tetra
 */
PRECICE_API void precicec_setMeshTetrahedron(
    const char *meshName,
    int         firstVertexID,
    int         secondVertexID,
    int         thirdVertexID,
    int         fourthVertexID);

/**
 * @brief Sets multiple mesh tetrahedra from vertex IDs
 *
 * @param[in] meshName the name of the mesh to add the tetrahedra to
 * @param[in] size the amount of tetrahedra to set
 * @param[in] ids an array containing 4*size vertex IDs
 *
 * @pre vertices in ids were added to the mesh with the name \p meshName
 */
PRECICE_API void precicec_setMeshTetrahedra(
    const char *meshName,
    int         size,
    const int * ids);

///@}

///@name Data Access
///@{

/// @copydoc precice::Participant::requiresInitialData()
PRECICE_API int precicec_requiresInitialData();

/**
 * @brief Writes vector data values given as block.
 *
 * The block must contain the vector values in the following form:
 * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
 * the number of vector values. In 2D, the z-components are removed.
 *
 * @param[in] meshName the name of the mesh
 * @param[in] dataName the name of the data to be written.
 * @param[in] size Number of indices, and number of values * dimensions.
 * @param[in] values Values of the data to be written.
 *
 * @see Participant::writeData
 */
PRECICE_API void precicec_writeData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const int *   valueIndices,
    const double *values);

/**
 * @brief Reads vector data values given as block.
 *
 * The block contains the vector values in the following form:
 * values = (d0x, d0y, d0z, d1x, d1y, d1z, ...., dnx, dny, dnz), where n is
 * the number of vector values. In 2D, the z-components are removed.
 *
 * @param[in] meshName the name of the mesh
 * @param[in] dataName the name of the data to be read.
 * @param[in] size  Number of indices, and number of values * dimensions.
 * @param[in] valueIndices Indices (from setReadPosition()) of data values.
 * @param[in] relativeReadTime Point in time where data is read relative to the beginning of the current time step.
 * @param[in] values Values of the data to be read.
 *
 * @see Participant::readData
 */
PRECICE_API void precicec_readData(
    const char *meshName,
    const char *dataName,
    int         size,
    const int * valueIndices,
    double      relativeReadTime,
    double *    values);

///@}

///@name Just-in-time mapping
///@{

/**
 * @brief Writes data using just-in-time data mapping
 * See @see precice::Participant::writeAndMapData() and the just-in-time mapping doxygen section for more information.
 *
 * @param[in] meshName Name of the mesh to write and map the data to. Typically a received mesh
 * @param[in] dataName Name of the data to field on this mesh.
 * @param[in] size Number of vertices to write
 * @param[in] coordinates Pointer to the coordinates where we write the data. Needs to have \p size x \ref precicec_getMeshDimensions( \p meshName ) entries
 * @param[in] values Pointer to the data values we want to pass to preCICE. Needs  to have \p size x \ref precicec_getDataDimensions( \p meshName , \p dataName ) entries
 *
 * @copydoc precice::Participant::writeAndMapData()
 */
PRECICE_API void precicec_writeAndMapData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const double *coordinates,
    const double *values);

/**
 * @brief Reads data using just-in-time data mapping.
 * See @see precice::Participant::mapAndReadData() and the just-in-time mapping doxygen section for more information.
 *
 * @param[in] meshName Name of the mesh to write and map the data to. Typically a received mesh
 * @param[in] dataName Name of the data to field on this mesh.
 * @param[in] size Number of vertices to write
 * @param [in] coordinates Pointer to the coordinates where we read the data. Needs to have \p size x \ref precicec_getMeshDimensions( \p meshName ) entries
 * @param[in] relativeReadTime Same as in @see precice::Participant::readData()
 * @param[in/out] values Pointer to the values to be filled by preCICE. Needs to have \p size x \ref precicec_getDataDimensions( \p meshName , \p dataName ) entries
 *
 * @copydoc precice::Participant::mapAndReadData()
 */
PRECICE_API void precicec_mapAndReadData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const double *coordinates,
    double        relativeReadTime,
    double *      values);
///@}

///@name Direct mesh access
///@{

/**
 * @brief See precice::Participant::setMeshAccessRegion().
 */
PRECICE_API void precicec_setMeshAccessRegion(
    const char *  meshName,
    const double *boundingBox);

/**
 * @brief See precice::Participant::getMeshVertexIDsAndCoordinates().
 */
PRECICE_API void precicec_getMeshVertexIDsAndCoordinates(
    const char *meshName,
    const int   size,
    int *       ids,
    double *    coordinates);

///@}

/** @name Experimental Gradient Data
 * These API functions are \b experimental and may change in future versions.
 */
///@{

/// @copydoc precice::Participant::requiresGradientDataFor
PRECICE_API int precicec_requiresGradientDataFor(const char *meshName,
                                                 const char *dataName);

/// @copydoc precice::Participant::writeGradientData
PRECICE_API void precicec_writeGradientData(
    const char *  meshName,
    const char *  dataName,
    int           size,
    const int *   valueIndices,
    const double *gradients);

///@}

/// @name User-defined profiling
///@{

/// @copydoc precice::Participant::startProfilingSection()
PRECICE_API void precicec_startProfilingSection(const char *sectionName);

/// @copydoc precice::Participant::stopLastProfilingSection()
PRECICE_API void precicec_stopLastProfilingSection();

///@}

/**
 * @brief Returns information on the version of preCICE.
 *
 * Returns a semicolon-separated C-string containing:
 *
 * 1) the version of preCICE
 * 2) the revision information of preCICE
 * 3) the configuration of preCICE including MPI, PETSC, PYTHON
 */
PRECICE_API const char *precicec_getVersionInformation();

#ifdef __cplusplus
}
#endif
