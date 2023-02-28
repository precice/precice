#pragma once

#include "precice/export.h"

/** @file
 * This file contains a Fortran 77 compatible interface written in C/C++.
 *
 *
 * Every method has a Fortran syntax equivalent in the method comment, and a
 * listing for input and output variables. A variable can be input and output
 * at the same time.
 */

#ifdef __cplusplus
extern "C" {
#endif

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
 * @copydoc precice::SolverInterface::SolverInterface()
 *
 */
PRECICE_API void precicef_create_(
    const char *participantName,
    const char *configFileName,
    const int * solverProcessIndex,
    const int * solverProcessSize,
    int         lengthAccessorName,
    int         lengthConfigFileName);

/**
 * Fortran syntax:
 * precicef_initialize( DOUBLE PRECISION timstepLengthLimit )
 *
 * IN:  -
 * OUT: timestepLengthLimit
 *
 * @copydoc precice::SolverInterface::initialize()
 *
 */
PRECICE_API void precicef_initialize_(double *timestepLengthLimit);

/**
 * Fortran syntax:
 * precicef_advance( DOUBLE PRECISION timstepLengthLimit )
 *
 * IN:  timestepLengthLimit
 * OUT: timestepLengthLimit
 *
 * @copydoc precice::SolverInterface::advance()
 *
 */
PRECICE_API void precicef_advance_(double *timestepLengthLimit);

/**
 * Fortran syntax:
 * precicef_finalize();
 *
 * @copydoc precice::SolverInterface::finalize()
 *
 */
PRECICE_API void precicef_finalize_();

/**
 * Fortran syntax:
 * precicef_get_dims( INTEGER dimensions )
 *
 * IN:  -
 * OUT: dimensions
 *
 * @copydoc precice::SolverInterface::getDimensions()
 *
 */
PRECICE_API void precicef_get_dims_(int *dimensions);

/**
 * Fortran syntax:
 * precicef_is_coupling_ongoing( INTEGER isOngoing )
 *
 * IN:  -
 * OUT: isOngoing(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::isCouplingOngoing()
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
 * @copydoc precice::SolverInterface::isTimeWindowComplete()
 *
 */
PRECICE_API void precicef_is_time_window_complete_(int *isComplete);

/**
 * Fortran syntax:
 * precicef_requires_initial_data_(
 *   INTEGER   isRequired )
 *
 * IN:  -
 * OUT: isRequired(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::requiresInitialData()
 */
PRECICE_API void precicef_requires_initial_data_(
    int *isRequired);

/**
 * Fortran syntax:
 * precicef_requires_reading_checkpoint_(
 *   INTEGER   isRequired )
 *
 * IN:  -
 * OUT: isRequired(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::requiresReadingCheckpoint()
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
 * @copydoc precice::SolverInterface::requiresWritingCheckpoint()
 */
PRECICE_API void precicef_requires_writing_checkpoint_(
    int *isRequired);

/**
 * Fortran syntax:
 * precicef_has_mesh(
 *   CHARACTER meshName(*),
 *   INTEGER   hasMesh )
 *
 * IN:  meshName
 * OUT: hasMesh(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::hasMesh()
 *
 */
PRECICE_API void precicef_has_mesh_(
    const char *meshName,
    int *       hasMesh,
    int         lengthMeshName);

/**
 * Fortran syntax:
 * precicef_has_data(
 *   CHARACTER dataName(*),
 *   INTEGER   meshID,
 *   INTEGER   hasData)
 *
 * IN:  dataName
 * IN:  meshID
 * OUT: hasData(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::hasData()
 *
 */
PRECICE_API void precicef_has_data_(
    const char *dataName,
    const int * meshID,
    int *       hasData,
    int         lengthDataName);

/**
 * Fortran syntax:
 * precicef_has_data(
 *   INTEGER   meshID
 *   INTEGER   required)
 *
 * IN:  meshID
 * OUT: required(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::requiresMeshConnectivityFor()
 */
PRECICE_API void precicef_requires_mesh_connectivity_for_(
    const int *meshID,
    int *      required);

/**
 * Fortran syntax:
 * precicef_set_vertex(
 *   INTEGER          meshID,
 *   DOUBLE PRECISION position(dim),
 *   INTEGER          vertexID )
 *
 * IN:  meshID, position
 * OUT: vertexID
 *
 * @copydoc precice::SolverInterface::setMeshVertex()
 *
 */
PRECICE_API void precicef_set_vertex_(
    const int *   meshID,
    const double *position,
    int *         vertexID);

/**
 * Fortran syntax:
 * precicef_get_mesh_vertex_size(
 *   INTEGER meshID,
 *   INTEGER meshSize )
 *
 * IN:  meshID
 * OUT: meshSize
 *
 * @copydoc precice::SolverInterface::getMeshVertexSize()
 *
 */
PRECICE_API void precicef_get_mesh_vertex_size_(
    const int *meshID,
    int *      meshSize);

/**
 * Fortran syntax:
 * precicef_set_vertices(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   DOUBLE PRECISION positions(dim*size),
 *   INTEGER          positionIDs(size) )
 *
 * IN:  meshID, size, positions
 * OUT: positionIDs
 *
 * @copydoc precice::SolverInterface::setMeshVertices()
 *
 */
PRECICE_API void precicef_set_vertices_(
    const int *meshID,
    const int *size,
    double *   positions,
    int *      positionIDs);

/**
 * Fortran syntax:
 * precicef_set_edge(
 *   INTEGER meshID,
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID )
 *
 * IN:  meshID, firstVertexID, secondVertexID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshEdge()
 *
 */
PRECICE_API void precicef_set_edge_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID);

/**
 * Fortran syntax:
 * precicef_set_mesh_edges_(
 *   INTEGER meshID,
 *   INTEGER size,
 *   INTEGER vertices(size*2) )
 *
 * IN:  meshID, size, vertices
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshEdges()
 *
 */
PRECICE_API void precicef_set_mesh_edges_(
    const int *meshID,
    const int *size,
    const int *vertices);

/**
 * Fortran syntax:
 * precicef_set_triangle_(
 *   INTEGER meshID,
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID )
 *
 * IN:  meshID, firstVertexID, secondVertexID, thirdVertexID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTriangle()
 *
 */
PRECICE_API void precicef_set_triangle_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID);

/**
 * Fortran syntax:
 * precicef_set_mesh_triangles_(
 *   INTEGER meshID,
 *   INTEGER size,
 *   INTEGER vertices(size*3) )
 *
 * IN:  meshID, size, vertices
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTriangles()
 *
 */
PRECICE_API void precicef_set_mesh_edges_(
    const int *meshID,
    const int *size,
    const int *vertices);

/**
 * Fortran syntax:
 * precicef_set_quad_(
 *   INTEGER meshID,
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID,
 *   INTEGER fourthVertexID )
 *
 * IN:  meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshQuad()
 *
 */
PRECICE_API void precicef_set_quad_(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID,
    const int *fourthVertexID);

/**
 * Fortran syntax:
 * precicef_set_mesh_quads(
 *   INTEGER meshID,
 *   INTEGER size,
 *   INTEGER vertices(size*4) )
 *
 * IN:  meshID, size, vertices
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshQuads()
 *
 */
PRECICE_API void precicef_set_mesh_quads_(
    const int *meshID,
    const int *size,
    const int *vertices);

/**
 * Fortran syntax:
 * precicef_set_tetrahedron(
 *   INTEGER meshID,
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID,
 *   INTEGER fourthVertexID )
 *
 * IN:  meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTetrahedron()
 *
 */
PRECICE_API void precicef_set_tetrahedron(
    const int *meshID,
    const int *firstVertexID,
    const int *secondVertexID,
    const int *thirdVertexID,
    const int *fourthVertexID);

/**
 * Fortran syntax:
 * precicef_set_mesh_tetrahedra_(
 *   INTEGER meshID,
 *   INTEGER size,
 *   INTEGER vertices(size*4) )
 *
 * IN:  meshID, size, vertices
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTetrahedra()
 *
 */
PRECICE_API void precicef_set_mesh_tetrahedra_(
    const int *meshID,
    const int *size,
    const int *vertices);

/**
 * Fortran syntax:
 * precicef_write_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  dataID, size, valueIndices, values
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockVectorData
 *
 */
PRECICE_API void precicef_write_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * Fortran syntax:
 * precicef_write_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim) )
 *
 * IN:  dataID, valueIndex, dataValue
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeVectorData
 *
 */
PRECICE_API void precicef_write_vdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue);

/**
 * Fortran syntax:
 * precicef_write_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size) )
 *
 * IN:  dataID, size, valueIndices, values
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockScalarData
 *
 */
PRECICE_API void precicef_write_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * Fortran syntax:
 * precicef_write_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue )
 *
 * IN:  dataID, valueIndex, dataValue
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeScalarData
 *
 */
PRECICE_API void precicef_write_sdata_(
    const int *   dataID,
    const int *   valueIndex,
    const double *dataValue);

/**
 * Fortran syntax:
 * precicef_read_bvdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  dataID, size, valueIndices
 * OUT: values
 *
 * @copydoc precice::SolverInterface::readBlockVectorData
 *
 */
PRECICE_API void precicef_read_bvdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * Fortran syntax:
 * precicef_read_vdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim) )
 *
 * IN:  dataID, valueIndex
 * OUT: dataValue
 *
 * @copydoc precice::SolverInterface::readVectorData
 *
 */
PRECICE_API void precicef_read_vdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue);

/**
 * Fortran syntax:
 * precicef_read_bsdata(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size) )
 *
 * IN:  dataID, size, valueIndices
 * OUT: values
 *
 * @copydoc precice::SolverInterface::readBlockScalarData
 *
 */
PRECICE_API void precicef_read_bsdata_(
    const int *dataID,
    const int *size,
    int *      valueIndices,
    double *   values);

/**
 * Fortran syntax:
 * precicef_read_sdata(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue )
 *
 * IN:  dataID, valueIndex
 * OUT: dataValue
 *
 * @copydoc precice::SolverInterface::readScalarData
 *
 */
PRECICE_API void precicef_read_sdata_(
    const int *dataID,
    const int *valueIndex,
    double *   dataValue);

PRECICE_API void precicef_get_version_information_(
    char *versionInfo,
    int   lengthVersionInfo);

/** @name Experimental Data Access
 * These API functions are \b experimental and may change in future versions.
 */
///@{

/**
 * Fortran syntax:
 * precicef_requires_gradient_data_for_(
 *   INTEGER dataID,
 *   INTEGER required )
 *
 * IN:  dataID
 * OUT: required(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::requiresGradientDataFor
 */
PRECICE_API void precicef_requires_gradient_data_for_(const int *dataID, int *required);

/**
 * Fortran syntax:
 * precicef_write_sgradient_data_(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION gradientValues )
 *
 * IN:  dataID, valueIndex, gradientValues
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeScalarGradientData
 */
PRECICE_API void precicef_write_sgradient_data_(
    const int *   dataID,
    const int *   valueIndex,
    const double *gradientValues);

/**
 * Fortran syntax:
 * precicef_write_bsgradient_data_(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION gradientValues )
 *
 * IN:  dataID, size, valueIndices, gradientValues
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockScalarGradientData
 */
PRECICE_API void precicef_write_bsgradient_data_(
    const int *   dataID,
    const int *   size,
    const int *   valueIndices,
    const double *gradientValues);

/**
 * Fortran syntax:
 * precicef_write_vgradient_data_(
 *   INTEGER dataID,
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION gradientValues )
 *
 * IN:  dataID, valueIndex, gradientValues
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeVectorGradientData
 */
PRECICE_API void precicef_write_vgradient_data_(
    const int *   dataID,
    const int *   valueIndex,
    const double *gradientValues);

/**
 * Fortran syntax:
 * precicef_write_bvgradient_data_(
 *   INTEGER dataID,
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION gradientValues )
 *
 * IN:  dataID, size, valueIndices, gradientValues
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockVectorGradientData
 */
PRECICE_API void precicef_write_bvgradient_data_(
    const int *   dataID,
    const int *   size,
    const int *   valueIndices,
    const double *gradientValues);

/**
 * Fortran syntax:
 * precicef_set_mesh_access_region_(
 *   INTEGER          meshID,
 *   DOUBLE PRECISION bounding_box(dim*2))
 *
 * IN:  meshID, bounding_box
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshAccessRegion()
 */
PRECICE_API void precicef_set_mesh_access_region_(
    const int     meshID,
    const double *boundingBox);

/**
 * Fortran syntax:
 * precicef_get_mesh_vertices_and_IDs_(
 *   INTEGER          meshID,
 *   INTEGER          size,
 *   INTEGER          ids(size),
 *   DOUBLE PRECISION coordinates(dim*size))
 *
 * IN:  meshID, size
 * OUT: ids, coordinates
 *
 * @copydoc precice::SolverInterface::getMeshVerticesAndIDs()
 */
PRECICE_API void precicef_get_mesh_vertices_and_IDs_(
    const int meshID,
    const int size,
    int *     ids,
    double *  coordinates);

///@}

#ifdef __cplusplus
}
#endif
