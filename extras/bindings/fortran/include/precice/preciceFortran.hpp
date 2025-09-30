#pragma once

#include <precice/export.h>

/** @file
 * This file contains a Fortran 77 compatible interface written in C/C++.
 * The Fortran bindings are thin wrappers around the C++ API.
 *
 * Every method has a Fortran syntax equivalent in the method comment, and a
 * listing for input and output variables. A variable can be input and output
 * at the same time.
 */

#ifdef __cplusplus
extern "C" {
#endif

///@name Construction and Configuration
///@{

/**
 * Fortran syntax:
 * precicef_create(
 *   CHARACTER participantName(*),
 *   CHARACTER configFileName(*),
 *   INTEGER   solverProcessIndex,
 *   INTEGER   solverProcessSize )
 *
 * IN:  participantName, configFileName, solverProcessIndex, solverProcessSize
 * OUT: -
 *
 * @copydoc precice::Participant::Participant()
 *
 */
PRECICE_API void precicef_create_(
    const char *participantName,
    const char *configFileName,
    const int  *solverProcessIndex,
    const int  *solverProcessSize,
    int         participantNameLength,
    int         configFileNameLength);

/**
 * Fortran syntax:
 * precicef_create_with_communicator_(
 *   CHARACTER participantName(*),
 *   CHARACTER configFileName(*),
 *   INTEGER   solverProcessIndex,
 *   INTEGER   solverProcessSize,
 *   INTEGER   communicator )
 *
 * IN:  participantName, configFileName, solverProcessIndex, solverProcessSize,
 *      communicator
 * OUT: -
 *
 * @copydoc precice::Participant::Participant()
 *
 */
PRECICE_API void precicef_create_with_communicator_(
    const char *participantName,
    const char *configFileName,
    const int  *solverProcessIndex,
    const int  *solverProcessSize,
    const int  *communicator,
    int         participantNameLength,
    int         configFileNameLength);

///@}

/// @name Steering Methods
///@{

/**
 * Fortran syntax:
 * precicef_initialize()
 *
 * IN:  -
 * OUT: -
 *
 * @copydoc precice::Participant::initialize()
 *
 */
PRECICE_API void precicef_initialize_();

/**
 * Fortran syntax:
 * precicef_advance( DOUBLE PRECISION timeStepSize )
 *
 * IN:  timeStepSize
 * OUT: -
 *
 * @copydoc precice::Participant::advance()
 *
 */
PRECICE_API void precicef_advance_(const double *timeStepSize);

/**
 * Fortran syntax:
 * precicef_finalize();
 *
 * @copydoc precice::Participant::finalize()
 *
 */
PRECICE_API void precicef_finalize_();

///@}

///@name Implicit coupling
///@{

/**
 * Fortran syntax:
 * precicef_requires_reading_checkpoint_(
 *   INTEGER   isRequired )
 *
 * IN:  -
 * OUT: isRequired(1:true, 0:false)
 *
 * @copydoc precice::Participant::requiresReadingCheckpoint()
 */
PRECICE_API void precicef_requires_reading_checkpoint_(
    int *isRequired);

/**
 * Fortran syntax:
 * precicef_requires_writing_checkpoint_(
 *   INTEGER   isRequired )
 *
 * IN:  -
 * OUT: isRequired(1:true, 0:false)
 *
 * @copydoc precice::Participant::requiresWritingCheckpoint()
 */
PRECICE_API void precicef_requires_writing_checkpoint_(
    int *isRequired);

///@}

///@name Status Queries
///@{

/**
 * Fortran syntax:
 * precicef_get_mesh_dimensions_(
 *   CHARACTER meshName(*),
 *   INTEGER   dimensions)
 *
 * IN:  mesh, meshNameLength
 * OUT: dimensions
 *
 * @copydoc precice::Participant::getMeshDimensions()
 *
 */
PRECICE_API void precicef_get_mesh_dimensions_(
    const char *meshName,
    int        *dimensions,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_get_data_dimensions_(
 *   CHARACTER meshName(*),
 *   CHARACTER dataName(*),
 *   INTEGER   dimensions)
 *
 * IN:  mesh, data, meshNameLength, dataNameLength
 * OUT: dimensions
 *
 * @copydoc precice::Participant::getDataDimensions()
 *
 */
PRECICE_API void precicef_get_data_dimensions_(
    const char *meshName,
    const char *dataName,
    int        *dimensions,
    int         meshNameLength,
    int         dataNameLength);

/**
 * Fortran syntax:
 * precicef_is_coupling_ongoing( INTEGER isOngoing )
 *
 * IN:  -
 * OUT: isOngoing(1:true, 0:false)
 *
 * @copydoc precice::Participant::isCouplingOngoing()
 *
 */
PRECICE_API void precicef_is_coupling_ongoing_(int *isOngoing);

/**
 * Fortran syntax:
 * precicef_is_time_window_complete( INTEGER isComplete );
 *
 * IN:  -
 * OUT: isComplete(1:true, 0:false)
 *
 * @copydoc precice::Participant::isTimeWindowComplete()
 *
 */
PRECICE_API void precicef_is_time_window_complete_(int *isComplete);

/**
 * Fortran syntax:
 * precicef_get_max_time_step_size( DOUBLE PRECISION maxTimeStepSize )
 *
 * IN:  -
 * OUT: maxTimeStepSize
 *
 * @copydoc precice::Participant::getMaxTimeStepSize()
 *
 */
PRECICE_API void precicef_get_max_time_step_size_(double *maxTimeStepSize);

///@}

///@name Mesh Access
///@anchor precice-mesh-access
///@{

/**
 * Fortran syntax:
 * precicef_precicef_requires_mesh_connectivity_for_(
 *   CHARACTER meshName(*),
 *   INTEGER   required)
 *
 * IN:  mesh, meshNameLength
 * OUT: required(1:true, 0:false)
 *
 * @copydoc precice::Participant::requiresMeshConnectivityFor()
 */
PRECICE_API void precicef_requires_mesh_connectivity_for_(
    const char *meshName,
    int        *required,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_reset_mesh_(
 *   CHARACTER meshName(*))
 *
 * IN:  mesh, meshNameLength
 *
 * @copydoc precice::Participant::resetMesh()
 */
PRECICE_API void precicef_reset_mesh_(
    const char *meshName,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_vertex(
 *   CHARACTER        meshName(*),
 *   DOUBLE PRECISION coordinates(dim),
 *   INTEGER          id )
 *
 * IN:  mesh, coordinates, meshNameLength
 * OUT: id
 *
 * @copydoc precice::Participant::setMeshVertex()
 *
 */
PRECICE_API void precicef_set_vertex_(
    const char   *meshName,
    const double *coordinates,
    int          *id,
    int           meshNameLength);

/**
 * Fortran syntax:
 * precicef_get_mesh_vertex_size(
 *   CHARACTER meshName(*),
 *   INTEGER meshSize )
 *
 * IN:  mesh, meshNameLength
 * OUT: meshSize
 *
 * @copydoc precice::Participant::getMeshVertexSize()
 *
 */
PRECICE_API void precicef_get_mesh_vertex_size_(
    const char *meshName,
    int        *meshSize,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_vertices(
 *   CHARACTER        meshName(*),
 *   INTEGER          size,
 *   DOUBLE PRECISION coordinates(dim*size),
 *   INTEGER          ids(size) )
 *
 * IN:  mesh, size, coordinates, meshNameLength
 * OUT: ids
 *
 * @copydoc precice::Participant::setMeshVertices()
 *
 */
PRECICE_API void precicef_set_vertices_(
    const char *meshName,
    const int  *size,
    double     *coordinates,
    int        *ids,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_edge(
 *   CHARACTER meshName(*),
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID )
 *
 * IN:  mesh, firstVertexID, secondVertexID, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::setMeshEdge()
 *
 */
PRECICE_API void precicef_set_edge_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_mesh_edges_(
 *   CHARACTER meshName(*),
 *   INTEGER size,
 *   INTEGER ids(size*2) )
 *
 * IN:  mesh, size, ids, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::setMeshEdges()
 *
 */
PRECICE_API void precicef_set_mesh_edges_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_triangle_(
 *   CHARACTER meshName(*),
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID )
 *
 * IN:  mesh, firstVertexID, secondVertexID, thirdVertexID, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::setMeshTriangle()
 *
 */
PRECICE_API void precicef_set_triangle_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_mesh_triangles_(
 *   CHARACTER meshName(*),
 *   INTEGER size,
 *   INTEGER ids(size*3) )
 *
 * IN:  mesh, size, ids, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::setMeshTriangles()
 *
 */
PRECICE_API void precicef_set_mesh_triangles_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_quad_(
 *   CHARACTER meshName(*),
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID,
 *   INTEGER fourthVertexID )
 *
 * IN:  mesh, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::setMeshQuad()
 *
 */
PRECICE_API void precicef_set_quad_(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    const int  *fourthVertexID,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_mesh_quads(
 *   CHARACTER meshName(*),
 *   INTEGER size,
 *   INTEGER ids(size*4) )
 *
 * IN:  mesh, size, ids, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::setMeshQuads()
 *
 */
PRECICE_API void precicef_set_mesh_quads_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_tetrahedron(
 *   CHARACTER meshName(*),
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID,
 *   INTEGER fourthVertexID )
 *
 * IN:  mesh, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::setMeshTetrahedron()
 *
 */
PRECICE_API void precicef_set_tetrahedron(
    const char *meshName,
    const int  *firstVertexID,
    const int  *secondVertexID,
    const int  *thirdVertexID,
    const int  *fourthVertexID,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_set_mesh_tetrahedra_(
 *   CHARACTER meshName(*),
 *   INTEGER size,
 *   INTEGER ids(size*4) )
 *
 * IN:  mesh, size, ids, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::setMeshTetrahedra()
 *
 */
PRECICE_API void precicef_set_mesh_tetrahedra_(
    const char *meshName,
    const int  *size,
    const int  *ids,
    int         meshNameLength);

///@}

///@name Data Access
///@{

/**
 * Fortran syntax:
 * precicef_requires_initial_data_(
 *   INTEGER   isRequired )
 *
 * IN:  -
 * OUT: isRequired(1:true, 0:false)
 *
 * @copydoc precice::Participant::requiresInitialData()
 */
PRECICE_API void precicef_requires_initial_data_(
    int *isRequired);

/**
 * Fortran syntax:
 * precicef_write_data(
 *   CHARACTER meshName(*),
 *   CHARACTER dataName(*),
 *   INTEGER size,
 *   INTEGER ids,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  mesh, data, size, ids, values, meshNameLength, dataNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::writeData
 *
 */
PRECICE_API void precicef_write_data_(
    const char *meshName,
    const char *dataName,
    const int  *size,
    int        *ids,
    double     *values,
    int         meshNameLength,
    int         dataNameLength);

/**
 * Fortran syntax:
 * precicef_read_data(
 *   CHARACTER meshName(*),
 *   CHARACTER dataName(*),
 *   INTEGER size,
 *   INTEGER ids,
 *   DOUBLE PRECISION relativeReadTime,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  mesh, data, size, ids, meshNameLength, dataNameLength
 * OUT: values
 *
 * @copydoc precice::Participant::readData
 *
 */
PRECICE_API void precicef_read_data_(
    const char   *meshName,
    const char   *dataName,
    const int    *size,
    int          *ids,
    const double *relativeReadTime,
    double       *values,
    int           meshNameLength,
    int           dataNameLength);

///@}

///@name Just-in-time mapping
///@{

/**
 * @brief Writes data using just-in-time data mapping
 * See @see precice::Participant::writeAndMapData() and the just-in-time mapping doxygen section for more information.
 *
 * Fortran syntax:
 * precicef_map_and_write_data_(
 *   CHARACTER meshName(*),
 *   CHARACTER dataName(*),
 *   INTEGER size,
 *   DOUBLE PRECISION coordinates(dataDim*size),
 *   DOUBLE PRECISION values(meshDim*size),
 *   INTEGER meshNameLength,
 *   INTEGER dataNameLength)
 *
 * IN:  mesh, data, size, coordinates, values, meshNameLength, dataNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::writeAndMapData()
 */
PRECICE_API void precicef_write_and_map_data_(
    const char *meshName,
    const char *dataName,
    const int  *size,
    double     *coordinates,
    double     *values,
    int         meshNameLength,
    int         dataNameLength);

/**
 * @brief Reads data using just-in-time data mapping.
 * See @see precice::Participant::mapAndReadData() and the just-in-time mapping doxygen section for more information.
 *
 * Fortran syntax:
 * precicef_map_and_read_data_(
 *   CHARACTER meshName(*),
 *   CHARACTER dataName(*),
 *   INTEGER size,
 *   DOUBLE PRECISION coordinates(dataDim*size),
 *   DOUBLE PRECISION relativeReadTime,
 *   DOUBLE PRECISION values(meshDim*size),
 *   INTEGER meshNameLength,
 *   INTEGER dataNameLength)
 *
 * IN:  mesh, data, size, coordinates, relativeReadTime, meshNameLength, dataNameLength
 * OUT: values
 *
 * @copydoc precice::Participant::mapAndReadData()
 *
 */
PRECICE_API void precicef_map_and_read_data_(
    const char   *meshName,
    const char   *dataName,
    const int    *size,
    double       *coordinates,
    const double *relativeReadTime,
    double       *values,
    int           meshNameLength,
    int           dataNameLength);
///@}

///@name Direct mesh access
///@{

/**
 * Fortran syntax:
 * precicef_set_mesh_access_region_(
 *   CHARACTER        meshName(*),
 *   DOUBLE PRECISION bounding_box(dim*2))
 *
 * IN:  mesh, bounding_box, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::setMeshAccessRegion()
 */
PRECICE_API void precicef_set_mesh_access_region_(
    const char   *meshName,
    const double *boundingBox,
    int           meshNameLength);

/**
 * Fortran syntax:
 * precicef_get_mesh_vertex_ids_and_coordinates_(
 *   CHARACTER        meshName(*),
 *   INTEGER          size,
 *   INTEGER          ids(size),
 *   DOUBLE PRECISION coordinates(dim*size))
 *
 * IN:  mesh, size, meshNameLength
 * OUT: ids, coordinates
 *
 * @copydoc precice::Participant::getMeshVertexIDsAndCoordinates()
 */
PRECICE_API void precicef_get_mesh_vertex_ids_and_coordinates_(
    const char *meshName,
    const int  *size,
    int        *ids,
    double     *coordinates,
    int         meshNameLength);

/**
 * Fortran syntax:
 * precicef_reset_mesh_access_region_(
 *   CHARACTER        meshName(*),
 *   INTEGER          meshNameLength)
 *
 * IN:  mesh, meshNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::resetMeshAccessRegion()
 */
PRECICE_API void precicef_reset_mesh_access_region_(
    const char *meshName,
    int         meshNameLength);
///@}

/** @name Experimental Gradient Data
 * These API functions are \b experimental and may change in future versions.
 */
///@{

/**
 * Fortran syntax:
 * precicef_requires_gradient_data_for_(
 *   CHARACTER meshName(*),
 *   CHARACTER dataName(*),
 *   INTEGER required )
 *
 * IN:  dataID
 * OUT: required(1:true, 0:false)
 *
 * @copydoc precice::Participant::requiresGradientDataFor
 */
PRECICE_API void precicef_requires_gradient_data_for_(
    const char *meshName,
    const char *dataName, int *required,
    int meshNameLength,
    int dataNameLength);

/**
 * Fortran syntax:
 * precicef_write_gradient_data_(
 *   CHARACTER meshName(*),
 *   CHARACTER dataName(*),
 *   INTEGER size,
 *   INTEGER ids,
 *   DOUBLE PRECISION gradients )
 *
 * IN:  mesh, data, size, ids, gradients, meshNameLength, dataNameLength
 * OUT: -
 *
 * @copydoc precice::Participant::writeGradientData
 */
PRECICE_API void precicef_write_gradient_data_(
    const char   *meshName,
    const char   *dataName,
    const int    *size,
    const int    *ids,
    const double *gradients,
    int           meshNameLength,
    int           dataNameLength);

///@}

/// @name User-defined profiling
///@{

/// @copydoc precice::Participant::startProfilingSection()
PRECICE_API void precicef_start_profiling_section_(const char *sectionName, int sectionNameLength);

/// @copydoc precice::Participant::stopLastProfilingSection()
PRECICE_API void precicef_stop_last_profiling_section_();

///@}

PRECICE_API void precicef_get_version_information_(
    char *versionInfo,
    int   lengthVersionInfo);

#ifdef __cplusplus
}
#endif
