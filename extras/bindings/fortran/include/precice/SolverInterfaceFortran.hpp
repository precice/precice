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
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER   hasData)
 *
 * IN:  mesh, data, lengthMesh, lengthData
 * OUT: hasData(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::hasData()
 *
 */
PRECICE_API void precicef_has_data_(
    const char *mesh,
    const char *data,
    int *       hasData,
    int         lengthMesh,
    int         lengthData);

/**
 * Fortran syntax:
 * precicef_precicef_requires_mesh_connectivity_for_(
 *   CHARACTER mesh(*),
 *   INTEGER   required)
 *
 * IN:  mesh, lengthMesh
 * OUT: required(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::requiresMeshConnectivityFor()
 */
PRECICE_API void precicef_requires_mesh_connectivity_for_(
    const char *mesh,
    int *       required,
    int         lengthMesh);

/**
 * Fortran syntax:
 * precicef_set_vertex(
 *   CHARACTER        mesh(*),
 *   DOUBLE PRECISION position(dim),
 *   INTEGER          vertexID )
 *
 * IN:  mesh, position, lengthMesh
 * OUT: vertexID
 *
 * @copydoc precice::SolverInterface::setMeshVertex()
 *
 */
PRECICE_API void precicef_set_vertex_(
    const char *  mesh,
    const double *position,
    int *vertexID int meshLenght);

/**
 * Fortran syntax:
 * precicef_get_mesh_vertex_size(
 *   CHARACTER mesh(*),
 *   INTEGER meshSize )
 *
 * IN:  mesh, meshLenght
 * OUT: meshSize
 *
 * @copydoc precice::SolverInterface::getMeshVertexSize()
 *
 */
PRECICE_API void precicef_get_mesh_vertex_size_(
    const char *mesh,
    int *       meshSize,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_set_vertices(
 *   CHARACTER        mesh(*),
 *   INTEGER          size,
 *   DOUBLE PRECISION positions(dim*size),
 *   INTEGER          positionIDs(size) )
 *
 * IN:  mesh, size, positions, meshLenght
 * OUT: positionIDs
 *
 * @copydoc precice::SolverInterface::setMeshVertices()
 *
 */
PRECICE_API void precicef_set_vertices_(
    const char *mesh,
    const int * size,
    double *    positions,
    int *       positionIDs,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_set_edge(
 *   CHARACTER mesh(*),
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID )
 *
 * IN:  mesh, firstVertexID, secondVertexID, meshLenght
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshEdge()
 *
 */
PRECICE_API void precicef_set_edge_(
    const char *mesh,
    const int * firstVertexID,
    const int * secondVertexID,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_set_mesh_edges_(
 *   CHARACTER mesh(*),
 *   INTEGER size,
 *   INTEGER vertices(size*2) )
 *
 * IN:  mesh, size, vertices, meshLenght
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshEdges()
 *
 */
PRECICE_API void precicef_set_mesh_edges_(
    const char *mesh,
    const int * size,
    const int * vertices,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_set_triangle_(
 *   CHARACTER mesh(*),
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID )
 *
 * IN:  mesh, firstVertexID, secondVertexID, thirdVertexID, meshLenght
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTriangle()
 *
 */
PRECICE_API void precicef_set_triangle_(
    const char *mesh,
    const int * firstVertexID,
    const int * secondVertexID,
    const int * thirdVertexID,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_set_mesh_triangles_(
 *   CHARACTER mesh(*),
 *   INTEGER size,
 *   INTEGER vertices(size*3) )
 *
 * IN:  mesh, size, vertices, meshLenght
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTriangles()
 *
 */
PRECICE_API void precicef_set_mesh_edges_(
    const char *mesh,
    const int * size,
    const int * vertices,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_set_quad_(
 *   CHARACTER mesh(*),
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID,
 *   INTEGER fourthVertexID )
 *
 * IN:  mesh, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID, meshLenght
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshQuad()
 *
 */
PRECICE_API void precicef_set_quad_(
    const char *mesh,
    const int * firstVertexID,
    const int * secondVertexID,
    const int * thirdVertexID,
    const int * fourthVertexID,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_set_mesh_quads(
 *   CHARACTER mesh(*),
 *   INTEGER size,
 *   INTEGER vertices(size*4) )
 *
 * IN:  mesh, size, vertices, meshLenght
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshQuads()
 *
 */
PRECICE_API void precicef_set_mesh_quads_(
    const char *mesh,
    const int * size,
    const int * vertices,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_set_tetrahedron(
 *   CHARACTER mesh(*),
 *   INTEGER firstVertexID,
 *   INTEGER secondVertexID,
 *   INTEGER thirdVertexID,
 *   INTEGER fourthVertexID )
 *
 * IN:  mesh, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID, meshLenght
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTetrahedron()
 *
 */
PRECICE_API void precicef_set_tetrahedron(
    const char *mesh,
    const int * firstVertexID,
    const int * secondVertexID,
    const int * thirdVertexID,
    const int * fourthVertexID,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_set_mesh_tetrahedra_(
 *   CHARACTER mesh(*),
 *   INTEGER size,
 *   INTEGER vertices(size*4) )
 *
 * IN:  mesh, size, vertices, meshLenght
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshTetrahedra()
 *
 */
PRECICE_API void precicef_set_mesh_tetrahedra_(
    const char *mesh,
    const int * size,
    const int * vertices,
    int         meshLenght);

/**
 * Fortran syntax:
 * precicef_write_bvdata(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  mesh, data, size, valueIndices, values, meshLength, dataLength
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockVectorData
 *
 */
PRECICE_API void precicef_write_bvdata_(
    const char *mesh,
    const char *data,
    const int * size,
    int *       valueIndices,
    double *    values,
    int         meshLenght,
    int         dataLength);

/**
 * Fortran syntax:
 * precicef_write_vdata(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim) )
 *
 * IN:  mesh, data, valueIndex, dataValue, meshLength, dataLength
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeVectorData
 *
 */
PRECICE_API void precicef_write_vdata_(
    const char *  mesh,
    const char *  data,
    const int *   valueIndex,
    const double *dataValue,
    int           meshLenght,
    int           dataLength);

/**
 * Fortran syntax:
 * precicef_write_bsdata(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size) )
 *
 * IN:  mesh, data, size, valueIndices, values, meshLength, dataLength
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockScalarData
 *
 */
PRECICE_API void precicef_write_bsdata_(
    const char *mesh,
    const char *data,
    const int * size,
    int *       valueIndices,
    double *    values,
    int         meshLenght,
    int         dataLength);

/**
 * Fortran syntax:
 * precicef_write_sdata(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue )
 *
 * IN:  mesh, data, valueIndex, dataValue, meshLength, dataLength
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeScalarData
 *
 */
PRECICE_API void precicef_write_sdata_(
    const char *  mesh,
    const char *  data,
    const int *   valueIndex,
    const double *dataValue,
    int           meshLenght,
    int           dataLength);

/**
 * Fortran syntax:
 * precicef_read_bvdata(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(dim*size) )
 *
 * IN:  mesh, data, size, valueIndices, meshLength, dataLength
 * OUT: values
 *
 * @copydoc precice::SolverInterface::readBlockVectorData
 *
 */
PRECICE_API void precicef_read_bvdata_(
    const char *mesh,
    const char *data,
    const int * size,
    int *       valueIndices,
    double *    values,
    int         meshLenght,
    int         dataLength);

/**
 * Fortran syntax:
 * precicef_read_vdata(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue(dim) )
 *
 * IN:  mesh, data, valueIndex, meshLength, dataLength
 * OUT: dataValue
 *
 * @copydoc precice::SolverInterface::readVectorData
 *
 */
PRECICE_API void precicef_read_vdata_(
    const char *mesh,
    const char *data,
    const int * valueIndex,
    double *    dataValue,
    int         meshLenght,
    int         dataLength);

/**
 * Fortran syntax:
 * precicef_read_bsdata(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION values(size) )
 *
 * IN:  mesh, data, size, valueIndices, meshLength, dataLength
 * OUT: values
 *
 * @copydoc precice::SolverInterface::readBlockScalarData
 *
 */
PRECICE_API void precicef_read_bsdata_(
    const char *mesh,
    const char *data,
    const int * size,
    int *       valueIndices,
    double *    values,
    int         meshLenght,
    int         dataLength);

/**
 * Fortran syntax:
 * precicef_read_sdata(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION dataValue )
 *
 * IN:  mesh, data, valueIndex, meshLength, dataLength
 * OUT: dataValue
 *
 * @copydoc precice::SolverInterface::readScalarData
 *
 */
PRECICE_API void precicef_read_sdata_(
    const char *mesh,
    const char *data,
    const int * valueIndex,
    double *    dataValue,
    int         meshLenght,
    int         dataLength);

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
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER required )
 *
 * IN:  dataID
 * OUT: required(1:true, 0:false)
 *
 * @copydoc precice::SolverInterface::requiresGradientDataFor
 */
PRECICE_API void precicef_requires_gradient_data_for_(const char *mesh,
                                                      const char *data, int *required,
                                                      int meshLenght,
                                                      int dataLength);

/**
 * Fortran syntax:
 * precicef_write_sgradient_data_(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION gradientValues )
 *
 * IN:  mesh, data, valueIndex, gradientValues, meshLength, dataLength
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeScalarGradientData
 */
PRECICE_API void precicef_write_sgradient_data_(
    const char *  mesh,
    const char *  data,
    const int *   valueIndex,
    const double *gradientValues,
    int           meshLenght,
    int           dataLength);

/**
 * Fortran syntax:
 * precicef_write_bsgradient_data_(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION gradientValues )
 *
 * IN:  mesh, data, size, valueIndices, gradientValues, meshLength, dataLength
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockScalarGradientData
 */
PRECICE_API void precicef_write_bsgradient_data_(
    const char *  mesh,
    const char *  data,
    const int *   size,
    const int *   valueIndices,
    const double *gradientValues,
    int           meshLenght,
    int           dataLength);

/**
 * Fortran syntax:
 * precicef_write_vgradient_data_(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER valueIndex,
 *   DOUBLE PRECISION gradientValues )
 *
 * IN:  mesh, data, valueIndex, gradientValues, meshLength, dataLength
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeVectorGradientData
 */
PRECICE_API void precicef_write_vgradient_data_(
    const char *  mesh,
    const char *  data,
    const int *   valueIndex,
    const double *gradientValues,
    int           meshLenght,
    int           dataLength);

/**
 * Fortran syntax:
 * precicef_write_bvgradient_data_(
 *   CHARACTER mesh(*),
 *   CHARACTER data(*),
 *   INTEGER size,
 *   INTEGER valueIndices,
 *   DOUBLE PRECISION gradientValues )
 *
 * IN:  mesh, data, size, valueIndices, gradientValues, meshLength, dataLength
 * OUT: -
 *
 * @copydoc precice::SolverInterface::writeBlockVectorGradientData
 */
PRECICE_API void precicef_write_bvgradient_data_(
    const char *  mesh,
    const char *  data,
    const int *   size,
    const int *   valueIndices,
    const double *gradientValues,
    int           meshLenght,
    int           dataLength);

/**
 * Fortran syntax:
 * precicef_set_mesh_access_region_(
 *   CHARACTER        mesh(*),
 *   DOUBLE PRECISION bounding_box(dim*2))
 *
 * IN:  mesh, bounding_box, meshLenght
 * OUT: -
 *
 * @copydoc precice::SolverInterface::setMeshAccessRegion()
 */
PRECICE_API void precicef_set_mesh_access_region_(
    const char *  mesh,
    const double *boundingBox,
    int           meshLenght);

/**
 * Fortran syntax:
 * precicef_get_mesh_vertices_and_IDs_(
 *   CHARACTER        mesh(*),
 *   INTEGER          size,
 *   INTEGER          ids(size),
 *   DOUBLE PRECISION coordinates(dim*size))
 *
 * IN:  mesh, size, lengthMesh
 * OUT: ids, coordinates
 *
 * @copydoc precice::SolverInterface::getMeshVerticesAndIDs()
 */
PRECICE_API void precicef_get_mesh_vertices_and_IDs_(
    const char *mesh,
    const int   size,
    int *       ids,
    double *    coordinates,
    int         meshLenght);

///@}

#ifdef __cplusplus
}
#endif
